/*
 *-----------------------------------------------------------------------------
 *
 * Object power spectral density estimation using a maximum likelihood approach
 * Based on routines originally written by Laurent Mugnier
 *
 * References :
 * J.M. Conan, L.M. Mugnier, T. Fusco, V. Michau & G. Rousset (1998) : "Myopic
 * deconvolution of adaptie optics images by use of object and point spread
 * function power spectra" Appl. Opt., vol 37, p. 4614
 *
 * A. Blanc, J. Idier,& L.M. Mugnier (2000) : "Novel estimator for the aberrations
 * of a space telescope by phase diversity", Proc. SPIE 4013, 728-736
 *
 *-----------------------------------------------------------------------------
 */

struct psd_struct
{
  pointer tfi2;
  pointer tfh2;
  pointer norm;
  long    nozero;
}

func psd_model(x,a,&grad,deriv=)
/* DOCUMENT psd_model(x,a,&grad,deriv=)
 *  
 * object psd simple model (a power law). 
 * psd= K / ((x/x0)^alpha + 1), where x is the radial frequency (in pixels).
 *
 */ 
{
  if (deriv == []) deriv=0;
  
  np=dimsof(x)(2);
  x0 = a(2)^2 + 0.1;
  alpha = a(3)^2;
  
  f = a(1) - log((exp(x)/x0)^alpha + 1.);
  
  if (deriv) {
    grad=array(0.,np,3);
    grad(,1) = 1.0;
    grad(,2) = 2.* a(2) * (alpha * exp(x)^alpha/x0^(alpha+1)) / ((exp(x)/x0)^alpha + 1.);
    grad(,3) = 2.* a(3) * (-(exp(x)/x0)^alpha * log(exp(x)/x0) / ((exp(x)/x0)^alpha + 1.));
  }

  return f;
}

func psd_fit(object,x0=,alpha=,disp=,verbose=,nozero=)
/* DOCUMENT psd_fit(object,x0=,alpha=,disp=,verbose=,nozero=)
 *  
 * fits the psd of object with a simple model.
 * uses lmfit (developped by Eric Thiebaut).
 *
 * KEYWORDS :
 * x0     : estimate of x0
 * alpha  : estimate of alpha
 * nozero : flag to include or not point [1,1]
 *
 */ 
{
  if (x0==[])    x0    = 1.;
  if (alpha==[]) alpha = 3.;

  np = dimsof(object)(2);

  if (disp==[]) disp=0;
  if (verbose==[]) verbose=0;
  if (nozero==[]) nozero=0;
  
  psd = (abs(fft(object,1))^2);
  psd_avg = circavg(psd);

  if (nozero) {
    x    = indgen(np/2-1);
    lpsd = log(psd_avg(2:np/2));
    w    = array(1.,np/2-1);
  } else {
    x    = indgen(np/2)-1.+1.e-3;
    lpsd = log(psd_avg(1:np/2));
    w    = array(1.,np/2);
  }
  
  a = [log(psd_avg(1)), x0, alpha]; 
  lnx = log(x);
  w /= sum(1./w);
  
  r= lmfit(psd_model,lnx,a,lpsd,w,deriv=1,stdev=1,monte_carlo=500,correl=1);

  k     = exp(a(1));
  x0    = a(2)^2+0.1; 
  alpha = a(3)^2;

  if (verbose)
    write,format="Estimated parameters : k = %.2f, x0 = %.2f, alpha = %.2f\n",k,x0,alpha;

  psd_fitted = double(k/(1.+(roll(dist(np))/x0)^alpha));
  
  if (disp) {
    window,1;
    fma;
    logxy,1,1;
    plg,circavg(psd),indgen(dimsof(circavg(psd))(2)),color="red",marks=0,width=4;
    plg,circavg(psd_fitted),indgen(dimsof(circavg(psd_fitted))(2)),
      color="green",marks=0,width=4;
    pltitle,"Fitted PSD (green)";
  }

  return psd_fitted;

}

func psd_error(param,&gradient,extra)
/* DOCUMENT psd_error(param,&gradient,extra)
 *  
 *   sum [ ln(|ft_psf|^2.psd_obj + psd_noise) +
 *              |ft_im-ft_psf.ft_mobj|^2/(|ft_psf|^2.psd_obj + psd_noise) ]
 *
 *
 */ 
{
  clip,param,1.e-16;
  
  tfi2   = *extra.tfi2;
  tfh2   = *extra.tfh2;
  norm   = *extra.norm;
  nozero = extra.nozero;
  
  n = (dimsof(tfi2))(2);

  x  = indgen(n)-1+1.e-3;

  param *= norm;

  sigma2 = param(1);
  k      = param(2);
  x0     = param(3);
  alpha  = param(4);

  psdo = k/((x/x0)^alpha+1);
  psdi = tfh2*psdo+sigma2;

  gradient=0.*param; 
  dpsdi = (1.)*(psdi-tfi2)/(psdi)^2;
  
  lnx     = x*0.;
  lnx(1)  = log(1.);
  lnx(2:) = log(x(2:n)/x0);
    
  if (nozero==1) dpsdi(1) = 0.;

  gradient(1) = sum(dpsdi);
  gradient(2) = sum(dpsdi*tfh2*psdo/k);
  gradient(3) = sum(dpsdi*tfh2*psdo^2/k*(alpha/x0)*(x/x0)^alpha);
  gradient(4) = -sum(dpsdi*tfh2*psdo^2/k*(x/x0)^alpha*lnx);

  gradient *= norm;

  if (nozero==1) {
    return sum((1.*(log(psdi)+tfi2/psdi))(2:n));
  } else return sum(1.*(log(psdi)+tfi2/psdi));
}


func yoda_psdobj(image,psf,&varnoise,&psd_noise,&psd_obj,mobj=,nozero=,disp=,
                 object=,nbiter=,threshold=,verbose=,noise_psf=)
/* DOCUMENT yoda_psdobj(image,psf,&varnoise,&psd_noise,&psd_obj,mobj=,nozero=,disp=,
 *                      object=,nbiter=,threshold=,verbose=,noise_psf=)
 *  
 *
 * Estimates the psd of the object from an image and a psf .
 * the object psd follows a simple model :
 * psd_obj = k / ((x/x0)^alpha + 1), where x is the radial frequency (in pixels).
 *
 * KEYWORDS :
 * varnoise  : Estimated noise variance
 * psd_noise : Estimated noise psd
 * psd_obj   : Estimated object psd
 * mobj      : Mean object (0 by default)
 * nozero    : flag : include or not point (1,1) in fit
 * disp      : disp flag
 * object    : actual object (if known)
 * nbiter    : max number of iterations 
 * threshold : convergence threshold 
 * verbose   : verbose flag
 * noise_psf : flag to force noise estimation and substraction on |FT(PSF)|^2
 *
 */ 
{
  if (nbiter==[])    nbiter    = 100;
  if (neval==[])     neval     = 2*nbiter;
  if (threshold==[]) threshold = 1.e-7;
  if (disp==[])      disp      = 0;
  if (nozero==[])    nozero    = 0;
  if (noise_psf==[]) noise_psf = 0;
   
  if (dimsof(image)(2) != sqrt(numberof(image)))
    error,"The image must be a square array";

  szh = dimsof(psf);
  if ((dimsof(psf)(2) != sqrt(numberof(psf))) || (numberof(psf) != numberof(image)))
    error,"The psf must be a square array of same size than image";

  if (mobj==[]) {
    mobj = 0.;
    tf_mobj = 0.;
  } else {
    if (numberof(mobj) == 1) {
      tf_mobj = 0.*image; 
      tf_mobj(np/2,np/2) = mobj; 
    } else tf_mobj= numberof(object) * roll(fft(roll(mobj),1));
  }

  tfh  = fft(roll(psf),1);
  tfh2 = abs(tfh)^2;
  
  np = long(sqrt(numberof(image)));
  
  if (noise_psf) {
    psd_noise_psf = (sum(roll(tfh2)(1,1:np))+sum(roll(tfh2)(np,1:np))+ \
                     sum(roll(tfh2)(2:np-1,1))+sum(roll(tfh2)(2:np-1,np)))/(4*np);
    tfh2-=psd_noise_psf;
    tfh2=clip(tfh2,1.e-7,);
  }
  
  tfh2 = tfh2 / max(tfh2);
  tfi2 = abs((fft(image,1)) - tfh * tf_mobj)^2;
  tfi2 = tfi2 / max(tfi2) * sum(image)^2;

  tfi2_avg = circavg(tfi2);
  tfh2_avg = circavg(tfh2);
  
  psd_noise = (sum(roll(tfi2)(1,1:np))+sum(roll(tfi2)(np,1:np))+ \
               sum(roll(tfi2)(2:np-1,1))+sum(roll(tfi2)(2:np-1,np)))/(4*np);
  
  if ((mobj!=0) || nozero) k = tfi2_avg(2);    
  else k = tfi2_avg(1);

  norm = double([psd_noise,k, 1., 1.]); 

  param    = [1.0, 1.0, 2.0, 2.0];

  extras        = psd_struct();
  extras.tfi2   = &tfi2_avg;
  extras.tfh2   = &tfh2_avg;
  extras.norm   = &norm;
  extras.nozero = nozero;

  nmin = numberof(param);
  if (is_void(ndir)) ndir = 5;
  method_name = swrite(format="Limited Memory BFGS (VMLM with NDIR=%d)",ndir);
  psd_ws = op_vmlmb_setup(nmin, ndir, fmin=0.0);
  
  step = 0.0;
  task = 1;
  eval = iter = 0;
  elapsed = array(double, 3);
  timer, elapsed;
  cpu_start = elapsed(1);
  stop_loop = 0;
  
  while (stop_loop == 0) {
    if (task == 1) {
      /* Evaluate function. */
      f = psd_error(param,g,extras);
      ++eval;
    }
    
    iter+=1;

    /* Check for convergence. */
    if (task != 1 || eval == 1) {
      timer, elapsed;
      cpu = elapsed(1) - cpu_start;
      gnorm = sqrt(sum(g*g));
      
      if (iter == nbiter) {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of iterations";
        stop_loop = 1;
      }
      
      if (eval == neval) {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of evaluations";
        stop_loop = 1;
      }
    }
    
    task = op_vmlmb_next(param, f, g, psd_ws);
    iter = (*psd_ws(2))(7);
    step = (*psd_ws(3))(22);

    if (task > 2) {
      write,"End of The Minimization process";
      write,"---------------------------------------------------------------";
      write,"I reached the convergence threshold";
      
      stop_loop = 1;
    }
  
    if (stop_loop) {
      newstring = swrite(format="Value of the criterion after %d iterations : %23.2e",
                         iter,psd_error(param,grad_fin,extras));
      write,"---------------------------------------------------------------";
      write,newstring;
    }
  }

  param *= norm;

  psd_noise = param(1);
  varnoise  = psd_noise/(np^2);
  k         = param(2);
  x0        = param(3);
  alpha     = param(4);

  if (verbose) {
    write,format="Estimated parameters : k = %.2f, x0 = %.2f, p = %.2f, sigma2 = %.2f\n",
      k,x0,alpha,varnoise;
  }
  
  psd_obj = k/(1.+(roll(dist(np))/x0)^alpha);
  
  if (disp) {
    s = numberof(tfi2_avg);
    window,2;
    fma;
    logxy,1,1;
    plg,circavg(tfi2),indgen(s),marks=0,width=4;
    plg,circavg(psd_obj)*tfh2_avg+psd_noise,indgen(s),color="red",marks=0,width=4;
    plg,indgen(s)*0.+psd_noise,indgen(s),color="cyan",marks=0,width=4;
    plg,circavg(psd_obj),indgen(s),color="green",marks=0,width=4;
    plt,"PSDi",0.62,0.88,tosys=0;
    plt,"PSDi est.",0.62,0.86,tosys=0,color="red";
    plt,"PSDo est.",0.62,0.84,tosys=0,color="green";
    plt,"PSDn",0.62,0.82,tosys=0,color="cyan";

    if (object!=[]) {
      real_psdo=abs(fft(object,1))^2;
      real_psdo=real_psdo/max(real_psdo)*sum(object)^2;
      fitted_psd=psd_fit(object,verbose=1,nozero=nozero,disp=0);
      plg,circavg(fitted_psd),indgen(s),color="magenta",marks=0,width=4;
      plg,circavg(real_psdo),indgen(s),color="blue",marks=0,width=4;
      plt,"PSDo fit.",0.62,0.80,tosys=0,color="magenta";
      plt,"PSDo real.",0.62,0.78,tosys=0,color="blue";
    }
    pltitle,"Object PSD estimate";
  }
 
  return [k,x0,alpha];
}
