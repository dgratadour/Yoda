/*
 * op_deconv.i
 *
 *	Regularized deconvolution in Yorick by various iterative optimization
 *	methods.
 */

require, "fft_utils.i";

func mks(n)
{
  a = array(long, n, n);
  a(1,1) = 3;
  a(0,0) = 3;
  if (n >= 3) a(n+2:n*(n-1)-1:n+1) = 2;
  a(2:n*(n-1):n+1) = 1;
  a(n+1:n*n-1:n+1) = 1;
  return 0.25*a;
}




func op_deconv(data, psf, x, weight=, rescale=, positive=,
               maxeval=, maxiter=, ndirs=, verb=, output=, save=,
               frtol=, fatol=, factr=, method=,
               mu=, wiener=, bootstrap=, normalize=,
               entropy=, beta=, prior=, tiny=,
               tikhonov=, roughness=)
/* DOCUMENT op_deconv(data, psf)
       -or- op_deconv(data, psf, x)


       FIXME:
       
     Returns the deconvolution  of DATA by PSF; X  is the optional starting
     solution.   The  solution  is   obtained  by  minimizing  the  penalty
     function:

       psi(X) = LKL(X) + MU*RGL(X)
       
     where RGL is a regularization term (see below) and the likelihood term
     is:
       
       LKL(X) = 0.5*sum(WEIGHT*(model(X) - DATA)^2)

     where WEIGHT is given by the value of optional keyword WEIGHT (default
     is WEIGHT=1.0), PRIOR is given  by the value of optional keyword PRIOR
     (default is PRIOR=0.0) and model(X) is the model of the data:

       model(X) = inverse_fft(fft(PSF)*fft(X)))


     If keyword SAVE is true, then if first least significant bit of SAVE
     is set (i.e. SAVE&1 == 1) the model is saved in external variable
     op_deconv_model and/or if second least significant bit of SAVE is set
     (i.e. SAVE&2 == 2) the penalty is saved in external variable
     op_deconv_penalty:
        op_deconv_penalty(1) = MU
        op_deconv_penalty(2) = LKL
        op_deconv_penalty(3) = RGL

     Keyword METHOD:
       METHOD = 0   - LBFGSB
              = 1   - VMLMB (no preconditioning)
              = 2   - VMLMB with preconditioning
              = 3   - conjugate gradient (Fletcher-Reeves)
              = 4   - conjugate gradient (Polak-Ribiére)
              = 5   - conjugate gradient (Polak-Ribiére with positive BETA)
     the default is METHOD=2
     
     If keyword BOOTSTRAP is true and X is not specified, then the starting
     solution is the solution of:

       X0 = arg min sum((model(X) - DATA)^2) + sum(PRIOR*abs(fft(X))^2)
          = arg min   sum(abs(fft(PSF)*fft(X) - fft(DATA))^2)
                    + N*sum(PRIOR*abs(fft(X))^2)
          = inverse_fft(conj(fft(PSF))*fft(DATA)/(abs(fft(PSF))^2 + N*PRIOR))
            
      with N=numberof(DATA).
            
      The definitions of the forward and backward FFT's are:
        fft(X) = fft(X, -1)
        inverse_fft(X) = (1.0/numberof(X))*double(fft(X, +1))

      With these definitions, Parseval's theorem reads:
        sum(abs(X)^2) = (1.0/numberof(X))*sum(abs(fft(X))^2)

      Keyword MU can be used to set the level of regularization (by default
      MU=1,  if some regularization  is used;  MU=0 otherwise).   There are
      several types of regularization:

       - Wiener regularization is used if WIENER keyword is set. The Wiener
         regularization term is:
 
           RGL_WIENER(X) = sum(WIENER*abs(fft(X))^2)
           
         The regularized unconstrained Wiener solution is given by:
         
           fft(X_WIENER) = conj(MTF)*fft(DATA)/(abs(MTF)^2 + MU*WIENER)

         If no  initial X  is specified, X_WIENER  is used as  the starting
         solution to find the constrained solution unless keyword BOOTSTRAP
         is explicitely set to 0.

       - Maximum entropy regularization is  used if keyword ENTROPY is true
         (non-nil and non-zero).  In this case, the regularization term is:

           RGL_ENTROPY(X) = BETA*KL(X, PRIOR) + (1 - BETA)*KL(PRIOR, X)

         where KL() is the Kullback-Leibler distance between the solution X
         and an a priori or default solution PRIOR.  KL is defined by:

           KL(X, P) = sum(X*log(X/P) + P - X)

         Keyword PRIOR  can be used to  specify the a priori,  by default a
         uniform prior is used.
         
         Keyword  BETA  can  be  used  to  adjust  the  definition  of  the
         regularization  term (by  default, BETA=0.5).   The regularization
         term reduces to the negative  Shannon entropy if BETA=1 and to the
         negative Burg entropy if BETA=0.

       - Tikhonov regularization is used
      

         This solution is obtained  value
         of the WIENER keyword is the Fourier weigths:
       
     
   SEE ALSO: deconv_lbfgsb. */
{
  require, "optim_pack.i";

  /* Setup for optional outputs. */
  extern op_deconv_model, op_deconv_penalty;
  if (! save) {
    save_model = save_penalty = 0n;
  } else {
    save_model   = (save&1 ? 1n : 0n);
    save_penalty = (save&2 ? 1n : 0n);
  }

  /* Initialize FFT stuff and compute MTF = FFT(PSF). */
  local __fft_setup, __fft_number;
  dims = dimsof(data);
  __fft_init, dims;
  fft = __fft;
  mtf = fft(psf, -1);
  conj_mtf = conj(mtf);
  fft_scale = 1.0/numberof(data);

  /* Type of regularization. */
  regul = ((is_void(wiener) ? 0 : 1)
           + (anyof(tikhonov)  ? 2 : 0)
           + (roughness ? 4 : 0)
           + (entropy   ? 8 : 0));
  if (is_void(mu)) mu = (regul ? 1.0 : 0.0);
  xmin = 0.0;
  if (regul == 0) {
    /* No regularization. */
    op = __op_deconv;
  } else if (regul == 1) {
    /* Wiener regularization. */
    op = __op_deconv_wiener;
    if (mu != 1.0) wiener *= mu;
  } else if (regul == 2) {
    /* Tikhonov regularization. */
    op = __op_deconv_tikhonov;
  } else if (regul == 4) {
    /* Roughness regularization. */
    op = (roughness==2 ? __op_deconv_relative_roughness
                       : __op_deconv_absolute_roughness);
  } else if (regul == 8) {
    /* Regularization by maximum entropy. */
    op = __op_deconv_entropy;
    positive = 1n; /* force positivity for maximum entropy method */
    if (is_void(beta)) beta = 0.5;
    else if (beta < 0.0 || beta > 1.0) error, "bad value for BETA";
    if (is_void(tiny)) tiny = 1e-30;
    else if (tiny <= 0.0) error, "TINY must be strictly positive";
    xavg = (normalize ? 1.0/numberof(data) : avg(data)/sum(psf));
    if (is_void(prior))
      prior = array(xavg, dimsof(data)); /* Uniform prior by default. */
    xmin = tiny*xavg;
  } else {
    error, "only one of WIENER, TIKHONOV, ROUGHNESS or ENTROPY is allowed";
  }
  
  /* Starting solution. */
  if (is_void(x)) {
    if (bootstrap) {
      x = fft(data, -1);
      if (is_void(wiener)) x /= mtf;
      else x *= conj_mtf/(abs2(mtf) + numberof(data)*wiener);
      x = fft_scale*double(fft(x, +1));
    } else {
      x = array(avg(data)/sum(psf), dims);
    }
  }

  /* Initialization of workspace according to optimization method
     and options. */
  if (is_void(frtol)) frtol = 1e-8;
  if (is_void(fatol)) fatol = 0.0;
  if (is_void(factr)) factr = 1e7;
  if (is_void(method)) method = 2; /* default is VMLMB with preconditioning */
  if (positive || method == 0) {
    /* Apply constraints and make XMIN an array to speedup things and
       because this is required by LBFGSB. */
    xmin = array(xmin, dims);
    if (positive) x = max(x, xmin);
  }
  precond = 0n;
  if (is_void(ndirs)) ndirs = 5;
  if (method == 0) {
    /* Use L-BFGS-B method. */
    if (! is_func(op_lbfgsb))
      error, "you must have L-BFGS-B compiled in Yorick";
    bnd = array((positive ? 1 : 0), dims);
    pgtol = 0.0;
    ws = op_lbfgsb_setup(ndirs, x, xmin, xmin /* unused upper bound */,
                         bnd, factr, pgtol);
    get_msg = op_lbfgsb_msg;
  } else if (method == 1 || method == 2) {
    /* Use variable metric. */
    fmin = 0.0; // FIXME:
    ws = op_vmlmb_setup(numberof(data), ndirs, fmin, fatol=fatol, frtol=frtol);
    get_msg = op_vmlmb_msg;
    if (method == 2) precond = 1n;
  } else {
    error, "bad METHOD";
  }
  job = 1;

  /* Do we need to explicitely apply constraints (not needed fo LBFGSB)? */
  apply_constraints = (positive && method != 0);
  
  /* Prepare for verbose output. */
  if (verb) {
    if (is_void(output)) {
      prefix = " ";
    } else {
      if (structof(output) == string) {
        output = open(output, "a");
      } else if (typeof(output) != "text_stream") {
        error, "bad value for keyword OUTPUT";
      }
      prefix = "#";
    }
    write, output, format="%s%s\n%s%s\n", prefix,
      "ITER  EVAL   FFT          PENALTY         LIKELIHOOD      REGUL.",
      prefix,
      "----  ----  -----  ---------------------  -----------  -----------";
  }

  /* Optimization loop. */
  iter = eval = 0;
  for (;;) {
    local lkl, lkl_h, rgl, rgl_h, grd, h, prior_over_x, sum_x; /* ouputs */
    
    if (job == 1) {
      /* Compute function and gradient. */

      /* Apply constraints. */
      if (apply_constraints) x = max(x, xmin);

      /* Compute model. */
      fft_x = fft(x, -1);
      model = fft_scale*double(fft(mtf*fft_x, +1));
      scale = (rescale ? op_deconv_scale() : 1.0);
      if (rescale) {
        if ((scale = op_deconv_scale(data, model, weight)) != 1.0) {
          model *= scale;
        }
      } else {
        scale = 1.0;
      }

      /* Compute penalty. */
      op, 1;
      err = lkl + mu*rgl;
      ++eval;
    }

    // FIXME: hack to save best solution found so far
    local best_err;
    if (eval == 1 || err < best_err) {
      if (save_model) eq_nocopy, op_deconv_model, model;
      best_err = err;
      best_x = (scale == 1.0 ? x : scale*x);
    }

    /* Check for convergence and, possibly, print iteration values. */
    if (job != 1 || eval == 1) {
      if (job == 2 || job == 3) ++iter;
      if ((stop = (job >= 3))) {
        msg = get_msg(ws);
      } else if ((stop = (! is_void(maxiter) && iter >= maxiter))) {
        msg = swrite(format="warning: too many iterations (%d)", iter);
      } else if ((stop = (! is_void(maxeval) && eval >= maxeval))) {
        msg = swrite(format="warning: too many function evaluations (%d)",
                     eval);
      }
      if (verb && (stop || iter%verb == 0)) {
        write, output, format=" %4d  %4d  %5d  %-22.14e %-11.3e %-11.3e\n",
          iter, eval, __fft_number, err, lkl, rgl;
      }
      if (stop) {
        /* Return current solution. */
        if (verb || job != 3) write, output, format=prefix+"%s\n", msg;
        if (save_penalty) op_deconv_penalty = [mu, lkl, rgl];
        return best_x;
      }
    }

    /* Call optimizer. */
    if (method == 0) {
      job = op_lbfgsb(x, err, grd, ws);
    } else {
      local h, active;
      if (job != 1 || eval == 1) {
        if (positive) active = (x > xmin) | (grd < 0.0);
        if (precond) op, 0;
      }
      //info, h;
      job = op_vmlmb_next(x, err, grd, ws, active, h);
    }
  }
#if 0
    if (0) {
      if (is_void(weight)) {
        alpha = sum(psf*psf);
      } else {
        alpha = fft_scale*double(fft(fft(weight, -1)*fft(psf*psf, -1), +1));
      }
    } else {
      if (is_void(weight)) {
        alpha = fft_scale*double(fft(abs2(mtf), +1));
      } else {
        alpha = fft_scale*double(fft(fft(weight, -1)*abs2(mtf), +1));
      }
    }
    if (numberof(alpha) > 1) {
      window,7;fma;pli,alpha;palette,"stern.gp";pause,1;
    }
    if (entropy) {
      tmp = x*x;
      tmp /= ((1.0 - beta)*prior + beta*x + alpha*tmp);
      next_stage, x, psi, grd, task, ws, bnd, tmp;
      tmp = [];
    } else {
      if (first_time) {
        tmp_i = 0;
        tmp_stride = 1;
        tmp_dims = dimsof(x);
        tmp_n = numberof(tmp_dims);
        tmp = array(0.0, tmp_dims);
        for (k=2;k<=tmp_n;++k) {
          tmp_i = tmp_i*tmp_stride + tmp_dims(k)/2;
          tmp_stride *= tmp_dims(k);
        }
        tmp(tmp_i) = mu;
        tmp = ((smooth(tmp) - tmp)(tmp_i))^2;
        write, tmp;
        alpha = 1.0/(alpha + tmp);
      }
      next_stage, x, psi, grd, task, ws, bnd, alpha;
    }
    first_time = 0n;
#endif
}



/* KL entropy:
   

rgl = (2.0*beta - 1.0)*(p - x) + ((beta - 1.0)*p + beta*x)*log(x/p);
grd = (1.0 - beta)*(1.0 - p/x) + beta*log(x/p);
hes = ((1.0 - beta)*(p/x) + beta)/x;

u = p/x;
 l = log(u);
 rgl = (2.0*beta - 1.0)*(p - x) + ((1.0 - beta)*p - beta*x)*l;
 grd = (1.0 - beta)*(1.0 - u) - beta*l;
 hes = ((1.0 - beta)*u + beta)/x;

 h = 1/diag(H)
   = 1/(lkl_h + mu*rgl_h)
   = x/(lkl_h*x + mu*((1.0 - beta)*u + beta));
*/

func __op_deconv_entropy(job)
{
  extern h, lkl_h, rgl_h;
  extern lkl, rgl, grd; /* ouputs */
  extern prior_over_x;
  if (job == 1) {
    /* Compute likelihood term and its gradient. */
    local weight_r;
    r = model - data; /* anti-residuals */
    if (is_void(weight)) eq_nocopy, weight_r, r;
    else weight_r = weight*r;
    grd = double(fft(conj_mtf*fft((fft_scale*scale)*weight_r, -1), +1));
    lkl = 0.5*sum(weight_r*r);

    /* Compute regularization. */
    prior_over_x = prior/x;
    log_prior_over_x = log(prior_over_x);
    rgl = sum((2.0*beta - 1.0)*(prior - x) +
              ((1.0 - beta)*prior - beta*x)*log_prior_over_x);
    grd += mu*((1.0 - beta)*(1.0 - prior_over_x) - beta*log_prior_over_x);
  } else {
    /* Compute preconditioning. */
    if (is_void(lkl_h)) {
      lkl_h = (is_void(weigth) ? sum(psf*psf) :
               fft_scale*double(fft(fft(weigth, -1)*fft(psf*psf, -1), +1)));
    }
    h = x/(lkl_h*x + mu*((1.0 - beta)*prior_over_x + beta));
  }
}


func op_deconv_scale(data, model, weight)
/* DOCUMENT op_deconv_scale(data, model)
       -or- op_deconv_scale(data, model, weight)
     Computes and returns the scale ALPHA such that the following quantity
     is minimized:
        sum(WEIGHT*(ALPHA*MODEL - DATA)^2)
     or sum((ALPHA*MODEL - DATA)^2)          if WEIGHT is missing.
     
   SEE ALSO: op_deconv. */
{
  if (is_void(weight)) {
    if ((tmp = sum(model*model)) > 0.0) return sum(model*data)/tmp;
  } else if ((tmp = sum((weight_model = weight*model)*model)) > 0.0) {
    return sum(weight_model*data)/tmp;
  }
  return 1.0;
}

func __op_deconv
{
  extern lkl, rgl, grd; /* ouputs */
  
  /* Compute likelihood term and (FFT of) gradient. */
  local weight_r;
  r = model - data; /* anti-residuals */
  if (is_void(weight)) eq_nocopy, weight_r, r;
  else weight_r = weight*r;
  grd = double(fft(conj_mtf*fft((fft_scale*scale)*weight_r, -1), +1));
  lkl = 0.5*sum(weight_r*r);
  rgl = 0.0;
}

func __op_deconv_tikhonov
{
  extern lkl, rgl, grd, sum_x; /* ouputs */
  
  /* Compute likelihood term and (FFT of) gradient. */
  local weight_r;
  r = model - data; /* anti-residuals */
  if (is_void(weight)) eq_nocopy, weight_r, r;
  else weight_r = weight*r;
  grd = double(fft(conj_mtf*fft((fft_scale*scale)*weight_r, -1), +1));
  rgl = 0.5*sum(weight_r*r);
  weight_r = r = [];

  if (normalize) {
    if ((sum_x = sum(x)) > 0.0) {
      u = x*(1.0/sum_x);
      v = tikhonov*u;
      uv = sum(u*v);
      grd += (mu/sum_x)*(v - uv);
      rgl = 0.5*uv;
      u = v = [];
    } else {
      rgl = 0.0;
    }
  } else {
    v = tikhonov*x;
    rgl = 0.5*sum(v*x);
    grd += mu*v;
  }
}

func __op_deconv_roughness
{
  extern h, lkl, rgl, grd, sum_x; /* ouputs */
  
  /* Compute likelihood term and (FFT of) gradient. */
  local weight_r;
  r = model - data; /* anti-residuals */
  if (is_void(weight)) eq_nocopy, weight_r, r;
  else weight_r = weight*r;
  grd = double(fft(conj_mtf*fft((fft_scale*scale)*weight_r, -1), +1));
  rgl = 0.5*sum(weight_r*r);
  weight_r = r = [];

  if (precond && eval == 1) {
    /* Computes the diagonal of the Hessian. */
    dims = dimsof(x);
    dims(2:) = 5;
    r = array(double, dims);
    r((numberof(r) + 1)/2) = 1.0;
    r = smooth(r) - r;
    rgl_h = sum(r*r);
    if (is_void(weigth)) {
      lkl_h = sum(psf*psf);
    } else {
      lkl_h = fft_scale*double(fft(fft(weigth, -1)*fft(psf*psf, -1), +1));
    }
    h = 1.0/(lkl_h + mu*rgl_h);
  }
  if (normalize) {
    if ((sum_x = sum(x)) > 0.0) {
      u = x*(1.0/sum_x);
      r = smooth(u) - u;
      rr = sum(r*r);
      grd += (mu/sum_x)*(smooth(r) - r - rr);
      rgl = 0.5*rr;
    } else {
      rgl = 0.0;
    }
  } else {
    r = smooth(x) - x;
    rgl = 0.5*sum(r*r);
    grd += mu*(smooth(r) - r);
  }
}


func op_deconv_lcurve(data, psf, x, weight=, rescale=, positive=,
                      maxeval=, maxiter=, ndirs=, verb=, output=, save=,
                      frtol=, fatol=, factr=, method=,
                      mu=, wiener=, bootstrap=, normalize=,
                      entropy=, beta=, prior=, tiny=,
                      tikhonov=, roughness=, precond=,
                      color=, win=, symbol=)
{
  /* sort regularization weight in descending order */
  mu = mu(sort(mu)(::-1));
  result = array(double, 3, numberof(mu));
  for (i=1 ; i<=numberof(mu) ; ++i) {
    local op_deconv_penalty;
    x = op_deconv(data, psf, x, weight=weight, rescale=rescale,
                  positive=positive, maxeval=maxeval, maxiter=maxiter,
                  ndirs=ndirs, verb=verb, output=output, save=2,
                  frtol=frtol, fatol=fatol, factr=factr, method=method,
                  mu=mu(i), wiener=wiener, bootstrap=bootstrap,
                  normalize=normalize, entropy=entropy, beta=beta,
                  prior=prior, tiny=tiny, tikhonov=tikhonov,
                  roughness=roughness);
    q = op_deconv_penalty;
    if (! is_void(win)) {
      window, win, wait=1;
      plp, q(3), q(2), color=color, symbol=symbol;
      pause, 1;
    }
    result(,i) = q;
  }
  return result;
}
