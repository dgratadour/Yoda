/*
    Yoda: YOrick Deconvolution Algorithm based on penalized maximum likelihood
    Copyright (C) <2011>  <Damien Gratadour>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Reference :
    L.M. Mugnier, T. Fusco & J.M. Conan (2004) : "MISTRAL: a myopic
    edge-preserving image restoration method, with application to astronomical
    adaptive-optics-corrected long-exposure images" OSAA, vol 21, p. 1841
 */

require,"imutil.i";
require,"/home/brujo/yorick/Yoda/trunk/yoda_utils.i";
require,"/home/brujo/yorick/Yoda/trunk/psd_fit.i";
require,"/home/brujo/yorick/Yoda/trunk/OptimPack-1.3.2/yorick/OptimPack1.i";
require,"/home/brujo/yorick/Yoda/trunk/OptimPack-1.3.2/yorick/lbfgs.i";

struct yoda_struct
{
  long    myopic;
  long    pos_obj;
  long    pos_psf;
  long    markov;
  long    nbiter;
  long    neval;
  long    verbose;
  string  regul;
  double  scale_psd;
  double  scale;
  double  delta;
  double  scale_psf;
  pointer image;
  pointer object;
  pointer psf;
  pointer psf_est;
  pointer ft_psf;
  pointer ft_image;
  pointer variance;
  pointer mobj;
  pointer ft_mobj;
  pointer psd;
  pointer psd_psf;
  pointer ft_mpsf;
  pointer prev_obj;
  pointer norm;
};



func yoda_fidelity(&grad_psf,&grad_obj,ft_object,image,ft_psf,variance)
/* DOCUMENT 
 * criterion=yoda_fidelity(&grad_psf,&grad_obj,ft_object,image,ft_psf,variance)
 *
 * sum (1/var x |(h * o) - i|^2)
 *
 */ 
{
  
  dim2 = numberof(image);
  np = long(sqrt(dim2));         //image must be square
  
  HO = 1./dim2*(double(fft((ft_psf*ft_object),-1)));
  
  HOminusI = HO - image;

  grad_psf = 1./dim2*double(fft(conj(ft_object)*fft(HOminusI/variance,1),-1));

  grad_obj = 1./dim2*double(fft(conj(ft_psf)*fft(HOminusI/variance,1),-1));

  return .5 * sum((HOminusI)^2/variance);
}

func yoda_fidelityf(&grad_psf,&grad_obj,ft_object,ft_image,ft_psf,variance)
/* DOCUMENT yoda_fidelityf(&grad_psf,&grad_obj,ft_object,ft_image,ft_psf,variance)
 *
 *  sum (1/var x |(h * o) - i|^2) in fourier space
 *
 */ 
{
  
  dim2 = numberof(image);
  np = long(sqrt(dim2));         //image must be square
  
  HO       = (ft_psf*ft_object);

  HOminusI = HO - ft_image;

  grad_psf = 1. / dim2 / variance * roll(fft(conj(ft_object) * HOminusI,-1).re);

  grad_obj = 1. / dim2 / variance * fft(conj(ft_psf) * HOminusI ,-1).re;

  return .5 / dim2 / variance * sum(abs(HOminusI)^2);
}

func yoda_markov(object,&grad_obj,scale=,delta=)
/* DOCUMENT yoda_markov(object,grad_object[,scale=,delta=])
 *  
 * delta^2 . sum [|dx|/delta -ln(1+|dx|/delta) + |dy|/delta -ln(1+|dy|/delta)]
 * where |dx| (m,n) = (o(m,n)-o(m-1,n))/scale is the x component of the
 * object gradient normalized by scale.
 *
 * KEYWORDS :
 * scale      : Scale factor applied to Do (default value is 1)
 * delta      : Threshold on the gradient for switching between linear and
 *              quadratic behavour. When delta tends to infinity, the criterion
 *              becomes purely quadratic.
 */ 
{

  if (scale==[]) scale = 1.;
  if (delta==[]) delta = 1.;

  dx = (object-roll(object,[1, 0])) / (delta * scale);
  dy = (object-roll(object,[0, 1])) / (delta * scale);

  crit = (delta^2)*sum(abs(dx)-log(1.+abs(dx))+ abs(dy) - log(1.+abs(dy)));

  dx /= (1. + abs(dx));
  dy /= (1. + abs(dy));

  grad_obj = (delta / scale)* (dx-roll(dx,[-1,0])+dy-roll(dy,[0,-1]));

  return crit;
}                                   

func yoda_l1l2(object,&grad_obj,scale=,delta=)
/* DOCUMENT yoda_l1l2(object,grad_object[,scale=,delta=])
 *  
 * delta^2 . sum [|dx|/delta -ln(1+|dx|/delta) + |dy|/delta -ln(1+|dy|/delta)]
 * where |dx| (m,n) = [o(m,n)-o(m-1,n)]/scale is the x component of the
 * object gradient normalized by scale.
 *
 * KEYWORDS :
 * scale      : Scale factor applied to Do (default value is 1)
 * delta      : Threshold on the gradient for switching between linear and
 *              quadratic behavour. When delta tends to infinity, the criterion
 *              becomes purely quadratic.
 *
 */ 
{

  if (!is_set(scale)) scale = 1.;
  if (!is_set(delta)) delta = 1.;

  dx = (object-roll(object,[1,0])) / (delta * scale);
  dy = (object-roll(object,[0,1])) / (delta * scale);
  
  r =sqrt(dx^2+dy^2);

  crit = (delta^2)*sum(r-log(1.+r));

  dx /= (1. + r);
  dy /= (1. + r);

  grad_obj = (delta / scale) * (dx-roll(dx,[-1,0])+dy-roll(dy,[0,-1]));

  return crit;
}


func yoda_gauss(&grad_obj,ft_obj,psd,ft_mobj)
/* DOCUMENT yoda_gauss(&grad_obj,ft_obj,psd,ft_mobj)
 *  
 * sum [|ft_object-ft_mobj|^2/psd_obj)]
 *
 *
 */ 
{
  if (ft_mobj==[]) ft_mobj = 0.;
  if (psd==[]) error,"PSD missing.";

  dim2=numberof(ft_object);

  grad_obj  = fft((ft_obj-ft_mobj)/psd,-1).re;
  
  return 0.5 * sum((abs(ft_obj-ft_mobj)^2)/psd);

}

func yoda_error(param,&gradient,extra)
/* DOCUMENT yoda_error(param,&gradient,extra)
 *  
 *
 */ 
{
  pos_obj  = extra.pos_obj;
  myopic   = extra.myopic;
  image    = *extra.image;
  variance = *extra.variance;
  norm     = *extra.norm;

  param   *= norm;
  
  sz=dimsof(param)(2);
  if (pos_obj) object = param(,1:sz)^2;
  else object = param(,1:sz);

  if (myopic) {
    pos_psf          = extra.pos_psf;
    if (pos_psf) psf = param(,sz+1:)^2;
    else psf = param(,sz+1:);

    psd_psf     = *extra.psd_psf;
    ft_mpsf     = *extra.ft_mpsf;
    scale_psf   = extra.scale_psf;
    
    norm_psf = sum(psf);
    psf /= norm_psf;
    ft_psf = fft(roll(psf),1);
    
    crit_array=array(0.,3); // data + regul_obj + regul_psf
  } else {
    ft_psf   = *extra.ft_psf;
    crit_array=array(0.,2); // data + regul
  }
  
  // the object is asumed to be a square array
  np = long(sqrt(numberof(object)));
  ft_object = fft(object,1);

  regul_type = extra.regul;
  
  if (regul_type == "l1-l2") {
    // we use the L1L2 regularization
    crit_array(1) = yoda_fidelity(grad_psf,grad_obj,ft_object,image,ft_psf,variance);
    
    scale     = extra.scale;
    delta     = extra.delta;
    markov    = extra.markov;
    
    if (markov!=0) {
      crit_array(2)=yoda_markov(object,grad_regul,scale=scale,delta=delta);
    } else {
      crit_array(2)=yoda_l1l2(object,grad_regul,scale=scale,delta=delta);
    }
    
    grad_obj += grad_regul; 

  } else {
    // else the gaussian prior is chosen
    ft_image = *yoda_data.ft_image;
    crit_array(1) = yoda_fidelityf(grad_psf,grad_obj,ft_object,ft_image,ft_psf,variance);
    psd       = *extra.psd;
    scale_psd = extra.scale_psd;
    ft_mobj   = *extra.ft_mobj;
    crit_array(2) = scale_psd * yoda_gauss(grad_regul,ft_object,psd,ft_mobj);
    grad_obj += scale_psd * grad_regul;
  }
    
  // if myopic is non null, we add the regularization term on the psf
  if (myopic) {
    crit_array(3) = scale_psf * yoda_gauss(grad_psf_regul,ft_psf,psd_psf,ft_mpsf);
    grad_psf += scale_psf * roll(grad_psf_regul);
  }

  gradient=param*0.;
  if (myopic) {
    if (pos_obj) gradient(,1:sz) = 2. * param(,1:sz) * grad_obj;
    else gradient(,1:sz) = grad_obj;
    if (pos_psf) gradient(,sz+1:) = 2. * param(,sz+1:) * (grad_psf - sum(psf*grad_psf))/ norm_psf;
    else gradient(,sz+1:) = (grad_psf - sum(psf*grad_psf))/ norm_psf;
  } else {
    if (pos_obj) gradient = 2. * param * grad_obj;
    else gradient = grad_obj;
  }
  gradient *= norm;

  if (extra.verbose) {
    write,"-----------------------------------------------------------\n";
    mystr1 = swrite(format = "Fidelity : %.2f / Regul obj : %.2f",crit_array(1),crit_array(2));
    if (myopic) mystr1 += swrite(format = " / Regul PSF : %.2f\n",crit_array(3));
    else mystr1 += "\n";
    write,mystr1;
    write,format = "Min/Max gradient obj : %.4f / %.4f\n",min(gradient(,1:sz)),max(gradient(,1:sz));
    if (myopic) write,format = "Min/Max gradient PSF : %.4f / %.4f\n",
                  min(gradient(,sz+1:)),max(gradient(,sz+1:));
    write,"-----------------------------------------------------------\n";
  }
  return sum(crit_array);
}


func yoda(image=,psf=,variance=,guess=,nbiter=,psd=,positivity=,mobj=,scale_psd=,
          scale=,delta=,markov=,disp=,threshold=,verbose=,myopic=,psd_psf=,mpsf=,
          scale_psf=,script=)
/* DOCUMENT func yoda(image=,psf=,variance=,guess=,nbiter=,psd=,positivity=,mobj=,scale_psd=,
 *                    scale=,delta=,markov=,disp=,threshold=,verbose=,myopic=,psd_psf=,mpsf=,
 *                    scale_psf=,script=)
 *
 *
 * The deconvolution can be regularized (gaussian or linear-quadratic
 * regularization) or not (the criterion contains only a term of fidelity to
 * the data). A positivity constraint can be added by reparametrization.
 * The deconvolution can be classical (known psf) or myopic (the psf is estimated
 * with the object).
 *
 * KEYWORDS :
 * image      : image 
 * psf        : psf (same size than image)
 * variance   : noise variance in the image (scalar or same size than image).
 * guess      : guess on the restored object
 * nbiter     : Number max of iterations
 * psd        : object psd (case of a quadratic regul)
 * positivity : flag for positivity 1 : o = x^2                
 * mobj       : mean object (quadratic regul). Default avg(image)
 * scale      : Scale factor applied to the object gradient for L1L2 regul.
 *            : Default value = 1. of the order of the standard deviation of the object gradient.
 * delta      : Threshold to switch between linear and quadratic regul
 *              When delta tends to infinity the regul is  quadratic
 * markov     : flag to select the Markov L1-L2 method
 * disp       : display flag
 * threshold  : threshold of the minimization process (default value =1.e-7  / 1.e-16)
 * verbose    : verbose flag
 * myopic     : myopic flag
 * psd_psf    : psf psd (myopic only)
 * mpsf       : mean psf (myopic only)
 * scale_psf  : scale of the psf regularization term (myopic only)
 *
 */ 
{
  extern yoda_data,yoda_iter;
  extern yoda_param,dispok,verboseok;
  extern stop_yoda_loop;
  extern elapsed,cpu_start,yoda_task,yoda_eval,yoda_step,yoda_ws,yoda_f,yoda_g;
  
  // Initialization of the verbose flag
  if (!is_set(verbose)) verbose=0;
  if (!is_set(disp)) disp=0;
  if (!is_set(script)) script=0;
  
  // Initialization of the myopic flag
  if (!is_set(myopic)) myopic=0;
  
  if (!is_set(pos_psf)) pos_psf=0;
  
  // ... testing the size of the image
  if (dimsof(image)(2) != dimsof(image)(3)) error,"The image must be a square array !";
  
  // ... testing the size of the psf
  if ((dimsof(psf)(2) != dimsof(psf)(3)) || (numberof(psf) != numberof(image)))
    error,"The PSF must be a square array of same size than image.";
  
  // ... testing the size and content of variance
  if (variance==[]) variance = 1.;
  else {
    if ((numberof(variance) != 1) & (numberof(variance) != numberof(image)))
      error,"variance must be a scalar or of same size than image.";
    else if (min(variance) <= 0.) error,"variance must be strictly positive";
  }
  
  np = long(sqrt(numberof(image)));

  psf = clip(psf,1.e-16,);
  psf /= sum(psf);

  // Create a guess if undifined or test its size
  if (myopic) {
    if (guess==[]) {
      guess=_(clip(image,avg(variance),),psf);
    }
    if (dimsof(guess)(3) != 2*dimsof(guess)(2))
      error,"guess must be npix x (2xnpix) for myopic deconv";
  } else {
    if (guess==[]) guess = clip(image,avg(variance),);
    if (dimsof(guess)(3) != dimsof(guess)(2)) error,"guess must be a square array";
  }

  ft_image = fft(image,1);
  ft_psf   = fft(roll(psf),1); // The PSF is supposed centered
                 
  if (numberof(nbiter) == 0) nbiter = 500; // 100 iterations by default
  neval = 5 * nbiter;
  if (mobj==[]) {
    if (verbose) write,"mobj is set to the average of the image";
    mobj = avg(guess(,1:np));
  }
  
  // mobj is needed for the gaussian prior :
  if (numberof(mobj) == 1) { // mobj scalar
    ft_mobj=0.*image;
    ft_mobj(1,1)=avg(guess(,1:np));
  } else ft_mobj =fft(mobj,1);
  
  if (myopic) {
    // mpsf is needed to estimate ft_mpsf
    if (mpsf!=[]) {
      if ((min(mpsf) < 0.)) error,"mpsf must be > 0";
      ft_mpsf = fft(roll(mpsf),1);
    } else {
      mpsf = avg(psf);
      ft_mpsf=0.*image;
      ft_mpsf(1,1)=mpsf;
    }
  }

  if (positivity) guess = sqrt(guess);
  
  norm     = guess*0. + 1.;
  if (myopic) {
    /*
    norm_obj = sum(guess(,1:np))/sqrt(np);
    norm(,1:np)   = norm_obj;
    guess(,1:np) /= norm_obj;
    */
    norm_obj = sum(guess(,1:np)/sqrt(variance));
    // normalization factor to get better convergence properties
    norm(,np+1:)   = 1./norm_obj;
    guess(,np+1:) *= norm_obj;
  }
  
  param = guess;
  
  if ((psd!=[]) && (is_set(delta)))
    error,"choose between Gaussian (keyword psd) and L1L2 (keyword delta) regularization";

  // Initialization of the convergence threshold
  if (numberof(threshold) == 0)
    threshold=1e-16;
  
  // Initialization of the yoda structure containing all the parameters
  if (delta==[]) delta=0;
  if (scale==[]) scale=0;
  if (markov==[]) markov=0;
  if (scale_psd==[]) scale_psd=1;
  
  yoda_data=yoda_struct();
  
  yoda_data.image          = &image;
  yoda_data.ft_image       = &(fft(image));
  yoda_data.psf            = &psf;
  yoda_data.ft_psf         = &ft_psf;
  yoda_data.variance       = &variance;
  yoda_data.mobj           = &mobj;
  yoda_data.ft_mobj        = &ft_mobj;
  yoda_data.psd            = &psd;
  yoda_data.pos_obj        = positivity;
  yoda_data.scale_psd      = scale_psd;
  yoda_data.scale          = scale;
  yoda_data.delta          = delta;
  yoda_data.markov         = markov;
  yoda_data.myopic         = myopic;
  yoda_data.nbiter         = nbiter;
  yoda_data.neval          = neval;
  if (delta != 0) yoda_data.regul = "l1-l2";
  if (myopic) {
    yoda_data.psd_psf    = &psd_psf;
    yoda_data.ft_mpsf    = &ft_mpsf;
    yoda_data.scale_psf  = scale_psf;
    yoda_data.pos_psf    = positivity;
  }
  yoda_data.norm           = &norm;
  yoda_data.verbose        = 0;
  
  // Some initializations for the minimization process
  yoda_param = double(param); //optimPack compat ... needs to be double
  nmin = numberof(yoda_param);
  min_method = 0;
  
  if (min_method == 0) {
    /* Variable metric. */
    if (is_void(ndir)) ndir = 5;
    method_name = swrite(format="Limited Memory BFGS (VMLM with NDIR=%d)",
                         ndir);
    yoda_ws = op_vmlmb_setup(nmin, ndir, fmin=threshold);
  } else if (min_method < 0) {
    if (is_void(ndir)) ndir = 5;
    method_name = swrite(format="Limited Memory BFGS (LBFGS with NDIR=%d)",
                         ndir);
    yoda_ws = op_lbfgs_setup(nmin, ndir);
  } else if (min_method >= 1 && min_method <= 15) {
    /* Conjugate gradient. */
    ndir = 2;
    method_name = swrite(format="Conjugate Gradient (%s)",
                         ["Fletcher-Reeves", "Polak-Ribiere",
                          "Polak-Ribiere with non-negative BETA"](min_method&3));
    yoda_ws = optim_cgmn_setup(min_method);
  }
  
  stop_yoda_loop = 0;
  dispok         = disp;
  verboseok      = verbose;
  
  yoda_step = 0.0;
  yoda_task = 1;
  yoda_eval = yoda_iter = 0;
  elapsed   = array(double, 3);
  timer, elapsed;
  cpu_start = elapsed(1);

  if (dispok) {
    winkill,0;
    window,0,dpi=100,style="nobox.gs",wait=1;
    animate,1;
  }
  
  if (script) {
    do {
      yoda_loop,1;
    } while (stop_yoda_loop == 0);
  } else yoda_loop;
  
}


func yoda_update_comments(iter, eval, cpu,f,gnorm,step)
{
  newstring = swrite(format=" %5d %5d %10.3f     %+-10.6e  %-9.1e  %-9.1e",iter, eval, cpu,f , gnorm, step);
  if (iter == 1) {
    write," ITER    EVAL     CPU [s]        FUNC       GNORM     STEPLEN  ";
    write,"------  ------  ----------   ----------  -----------    -------";
    write,newstring;
  } else write,format="\r %s \n",newstring;
}

   

func yoda_loop(one)
{
  extern yoda_data,yoda_iter;
  extern yoda_param,dispok,verboseok;
  extern stop_yoda_loop;
  extern elapsed,cpu_start,yoda_task,yoda_eval,yoda_step,yoda_ws,yoda_f,yoda_g;
     
  if (yoda_task == 1) {
    /* Evaluate function. */
    yoda_f = yoda_error(yoda_param,yoda_g,yoda_data);
    ++yoda_eval;
  }

  if (myopic) {
    sz = long(sqrt(numberof(*yoda_data.image)));
    yoda_data.object  = &(yoda_param(,1:sz)^2);
    yoda_data.psf_est = &(yoda_param(,sz+1:)^2);
  } else {
    yoda_data.object  = &(yoda_param^2);
  }

  if (yoda_iter == 0) yoda_data.prev_obj = &(*yoda_data.object);
  
  yoda_iter+=1;
  //error;
  /* Check for convergence. */
  if (yoda_task != 1 || yoda_eval == 1) {
    timer, elapsed;
    cpu = elapsed(1) - cpu_start;
    gnorm = sqrt(sum(yoda_g*yoda_g));
    
    if (verboseok) yoda_update_comments,yoda_iter, yoda_eval, cpu,yoda_f,gnorm,yoda_step;
    
    if (yoda_iter == yoda_data.nbiter) {
      if (verboseok) {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of iterations";
      }
      stop_yoda_loop = 1;
    }

    if (yoda_eval == yoda_data.neval) {
      if (verboseok) {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of evaluations";
      }
      stop_yoda_loop = 1;
    }
    
    if (dispok) yoda_update_disp,yoda_data;
  }
  
  min_method = 0;
  
  if (min_method == 0) {
      yoda_task = op_vmlmb_next(yoda_param, yoda_f, yoda_g, yoda_ws);
      yoda_iter = (*yoda_ws(2))(7);
      yoda_step = (*yoda_ws(3))(22);
  } else {
    // compat with optim not yet working stuffs !
    if (min_method < 0) {
      yoda_task = op_lbfgs_next(double(yoda_param), yoda_f, yoda_g, yoda_ws);
      if (yoda_task == 2 || yoda_task == 3) ++yoda_iter;
      yoda_step = -1.0;
    } else {
      optim_cgmn, yoda_param, yoda_f, yoda_g, yoda_task, yoda_ws;
      yoda_iter = optim_cgmn_iter(yoda_ws);
      yoda_step = optim_cgmn_step(yoda_ws);
    }
  }
  
  if (yoda_task > 2) {
    write,"End of The Minimization process";
    write,"---------------------------------------------------------------";
    write,"I reached the convergence threshold";
      
    stop_yoda_loop = 1;
  }
  
  if (stop_yoda_loop) {
    newstring = swrite(format="Value of the criterion after %d iterations : %23.2e",
                       yoda_iter,yoda_error(yoda_param,grad_fin,yoda_data));
    write,"---------------------------------------------------------------";
    write,newstring;
    animate,0;
    //set_idler;
    return;
  }
  
  if (!one) {
    set_idler,yoda_loop;
  }
}

func pause(void)
{
  extern stop_yoda_loop;

  stop_yoda_loop = 1;
  animate,0;
}
func yoda_update_disp(yoda_data)
{
  np = long(sqrt(numberof(*yoda_data.image)));
  imrec = array(0.,2*np,2*np);
  imrec(1:np,1:np)   = clip(*yoda_data.image,0.,);
  if (myopic) imrec(np+1:,1:np) = *yoda_data.psf_est/max(*yoda_data.psf_est)*max(*yoda_data.image);
  else imrec(np+1:,1:np) = *yoda_data.psf/max(*yoda_data.psf)*max(*yoda_data.image);
  imrec(1:np,np+1:)  = *yoda_data.object/max(*yoda_data.object)*max(*yoda_data.image);
  imrec(np+1:,np+1:) = abs(*yoda_data.object - *yoda_data.prev_obj)^2;
  yoda_data.prev_obj = &(*yoda_data.object);
  window,0;fma;
  pli,imrec;
  plt,"Object",0.37,0.64,tosys=0,color="red",height=8;
  plt,"!Dobj",0.405,0.64,tosys=0,color="red",height=8;
  plt,"Image",0.37,0.625,tosys=0,color="red",height=8;
  plt,"PSF",0.405,0.625,tosys=0,color="red",height=8;
  mystr = swrite(format="Iteration #: %d / Error = %.2f",yoda_iter,yoda_f);
  plt,mystr,0.47,0.37,tosys=0,color="red",height=8;
  redraw;
}

/*
grad_oc = image*0.;
ft_object = fft(object,1);
ref=yoda_fidelityf(grad_psf,grad_obj,ft_object,ft_image,ft_psf,variance);
for (i=1;i<=numberof(image);i++) {
  obj2 = object;
  obj2(i) += 1;
  ft_object = fft(obj2,1);
  res=yoda_fidelityf(grad_psf,grad_obj2,ft_object,ft_image,ft_psf,variance);
  grad_oc(i) = res-ref;
 }
avg(grad_oc/grad_obj)
plmk,grad_oc(*),grad_obj(*)

grad_oc = psf*0.;
ft_object = fft(object,1);
ft_psf = fft(roll(psf),1);
ref=yoda_fidelityf(grad_psf,grad_obj,ft_object,ft_image,ft_psf,variance,const_psf);
for (i=1;i<=numberof(image);i++) {
  psf2 = psf;
  psf2(i) += 1.e-7;
  norm_psf2 = sum(psf2);
  psf2 /= norm_psf2;
  ft_psf = fft(roll(psf2),1);
  res=yoda_fidelityf(grad_psf2,grad_obj2,ft_object,ft_image,ft_psf,variance);
  grad_oc(i) = (res-ref)*1.e7;
 }
avg(grad_oc/(grad_psf-const_psf))
plmk,grad_oc(*),grad_psf(*)


grad_oc = image*0.;
ft_object = fft(object,1);
ref=yoda_gauss(grad_regul,ft_object,psd,ft_mobj);
for (i=1;i<=numberof(image);i++) {
  obj2 = object;
  obj2(i) += 1;
  ft_object = fft(obj2,1);
  res=yoda_gauss(grad_regul2,ft_object,psd,ft_mobj);
  grad_oc(i) = res-ref;
 }
avg(grad_oc/grad_regul) 
plmk,grad_oc(*),grad_regul(*)


grad_oc = image*0.;
ft_psf = fft(roll(psf),1);
ref=yoda_gauss(grad_regul,ft_psf,psd_psf,ft_mpsf);
for (i=1;i<=numberof(image);i++) {
  psf2 = psf;
  psf2(i) += 1.e-10;
  norm_psf2 = sum(psf2);
  psf2 /= norm_psf2;
  ft_psf = fft(roll(psf2),1);
  res=yoda_gauss(grad_regul2,ft_psf,psd_psf,ft_mpsf);
  grad_oc(i) = (res-ref)*1.e10;
 }
avg(grad_oc/roll(grad_regul)) 
plmk,grad_oc(*),roll(grad_regul)(*)


grad_oc = yoda_param*0.;
param = yoda_param;
sz=dimsof(param)(2);
ref=yoda_error(param,gradient,yoda_data);
for (i=1;i<=numberof(param);i++) {
  param = yoda_param;
  if (float(i)/(sz*sz) > 1) param(i) += 1.e-10;
  else param(i) += 1.;
  res=yoda_error(param,gradient2,yoda_data);
  if (float(i)/(sz*sz) > 1) grad_oc(i) = (res-ref)*1.e10;
  else grad_oc(i) = res-ref;
}


*/

