require,"yao.i";
require,"yoda.i";

func test_yoda_gaussian(nlevel)
{
  object=10.*fits_read("yoda.fits");
  nobj = 128;
  
  /*
  obj_bin = object(cum,cum)(::2,::2)(dif,dif);
  object *= 0.;
  object(32:32+63,32:32+63) = obj_bin;
  */
  /*
  nobj = 256;
  obj_big = array(0.,nobj,nobj);
  obj_big(nobj/2-128/2:nobj/2+128/2-1,nobj/2-128/2:nobj/2+128/2-1) = object;  
  object = obj_big;
  */
  
  psf=telfto(2.2,0.5,8.,0.11,0.013,nobj,returnpsf=1);
  psf=psf/sum(psf);

  image=rfftconvol(object,psf);
  psf=roll(psf);

  varnoise = max(image/nlevel)^2;

  noise       = random_n([2,nobj,nobj])*sqrt(varnoise);
  image_noise = image+ noise;
  
  psf_noise   = (psf*sum(image)+random_n([2,nobj,nobj])*sqrt(varnoise))/sum(image);

  write,"Estimating the object spectrum from the image";
  
  psdo=abs(fft(object,1))^2;
  psdo=clip(psdo,1.e-6,);
  //psdo=roll(psdo);
  psdo=psdo/max(psdo)*(sum(image_noise))^2;
  
  params=yoda_psdobj(image_noise,psf,varn_est,psd_noise,psd_obj,disp=1,
                    verbose=0,nozero=0,threshold=1.e-6,noise_psf=1,object=object);

  ratio = psd_noise/psd_obj;
  
  //computing wiener estimate
  myguess = fft((fft(image_noise)*fft(psf)/(abs(fft(psf))^2+ratio)),-1).re / numberof(image);
  myguess = roll(myguess);
  myguess /= sum(myguess);
  myguess *= sum(object);

  alpha = (max(ratio)>1. ? 1./max(ratio) : 1.);

  write,"Starting deconvolution with yoda ...";

  write,"Deconvolution with a gaussian prior";
  
  yoda,image=image_noise,psf=psf,variance=varnoise,guess=clip(image,1.e-6,),positivity=1,
    scale_psd = alpha ,psd=psdo,mobj=object,disp=1,script=1;
  
  known_psd = *yoda_data.object;
  known_psd = known_psd / avg(known_psd) * avg(object);

  psd_obj = psd_obj / max(psd_obj) * sum(image)^2;

  yoda,image=image_noise,psf=psf,variance=varnoise,guess=clip(image,1.e-6,),positivity=1,
    scale_psd= alpha ,psd=psd_obj,disp=1,script=1;
  
  est_psd = *yoda_data.object;
  est_psd = est_psd / avg(est_psd) * avg(object);

  winkill,5;
  window,5,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = est_psd-object;
  imrec(1:nobj,nobj+1:)  = known_psd-object;
  imrec(nobj+1:,nobj+1:) = myguess-object;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Error (Gaussian prior vs Wiener)";

  winkill,4;
  window,4,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = image_noise;
  imrec(1:nobj,nobj+1:)  = est_psd;
  imrec(nobj+1:,nobj+1:) = myguess;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Gaussian prior vs Wiener";
  
  write,"wiener error       :",sum((myguess-object)^2)/sum(object^2);
  write,"estimated psd error:",sum((est_psd-object)^2)/sum(object^2);
  write,"real psd error     :",sum((known_psd-object)^2)/sum(object^2);

  return imrec;
}

func test_yoda_l1l2(scale,delta)
{
  if (scale == []) scale = 100.;
  if (delta == []) delta = 1.;
  
  object=10.*fits_read("yoda.fits");
  /*
  obj_bin = object(cum,cum)(::2,::2)(dif,dif);
  object *= 0.;
  object(32:32+63,32:32+63) = obj_bin;
  */
  nobj = 256;
  obj_big = array(0.,nobj,nobj);
  obj_big(nobj/2-128/2:nobj/2+128/2-1,nobj/2-128/2:nobj/2+128/2-1) = object;  
  object = obj_big;
  
  psf=telfto(2.2,0.5,8.,0.11,0.013,nobj,returnpsf=1);
  psf=psf/sum(psf);
  psdo=abs(fft(object,1))^2;
  psdo=clip(psdo,max(psdo)*1.e-7,);
  psdo=roll(psdo);///max(psdo)*sum(object));

  image=rfftconvol(object,psf);
  psf=roll(psf);

  varnoise = max(image/50.)^2;

  image_noise=image+(random_n([2,nobj,nobj])-0.5)*sqrt(varnoise);
  psf_noise=(psf*sum(image)+random_n([2,nobj,nobj])*sqrt(varnoise))/sum(image);

  //computing wiener estimate
  myguess = fft((fft(image)*fft(psf)/(abs(fft(psf))^2+1./varnoise)),-1).re /
    numberof(image);
  myguess = roll(myguess);
  
  write,"Starting deconvolution with yoda ...";
  /*
  write,"Deconvolution with a L1L2 prior no positivity";
  yoda,image=image_noise,psf=psf,variance=varnoise,guess=myguess,
    positivity=0,scale=scale,delta=delta,disp=1,script=1;
  nopos = *yoda_data.object;
  nopos = nopos / avg(nopos) * avg(object);
  */
  
  write,"Deconvolution with a L1L2 prior + positivity";
  yoda,image=image_noise,psf=psf,variance=varnoise,guess=clip(image,1.e-7,),
    positivity=1,scale=scale,delta=delta,disp=1,script=1;
  wpos = *yoda_data.object;
  wpos = wpos / avg(wpos) * avg(object);
  /*
  write,"Deconvolution with a Markov L1L2 prior + positivity";
  yoda,image=image,psf=psf,variance=varnoise,guess=clip(image,1.e-7,),
    positivity=2,markov=1,scale=scale,delta=delta,disp=1, script=1;
  mrkv = *yoda_data.object;
  mrkv = mrkv / avg(mrkv) * avg(object);
  */
  mrkv = 0.;
  nopos = 0.;
  winkill,5;
  window,5,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = myguess-object;
  imrec(1:nobj,nobj+1:)  = wpos-object;
  imrec(nobj+1:,nobj+1:) = mrkv-object;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Error (L1-L2 prior)";

  winkill,4;
  window,4,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = myguess;
  imrec(1:nobj,nobj+1:)  = wpos;
  imrec(nobj+1:,nobj+1:) = mrkv;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Deconvolution with L1-L2 prior";
  
  write,"no positivity error:",sum((nopos-object)^2)/sum(object^2);
  write,"wiener error       :",sum((myguess-object)^2)/sum(object^2);
  write,"w/ positivity error:",sum((wpos-object)^2)/sum(object^2);
  write,"markov error       :",sum((mrkv-object)^2)/sum(object^2);
  
  return imrec;

}


func test_yoda_gaussian_myopic(nlevel)
{
  object=10.*fits_read("yoda.fits");
  nobj = 128;
  
  /*
  obj_bin = object(cum,cum)(::2,::2)(dif,dif);
  object *= 0.;
  object(32:32+63,32:32+63) = obj_bin;
  */

  /*
  nobj = 256;
  obj_big = array(0.,nobj,nobj);
  obj_big(nobj/2-128/2:nobj/2+128/2-1,nobj/2-128/2:nobj/2+128/2-1) = object;  
  object = obj_big;
  */
  
  psf=telfto(2.2,0.5,8.,0.11,0.013,nobj,returnpsf=1);
  psf=psf/sum(psf);

  nmodes  = 100;
  
  cent    = nobj/2+0.5;
  pupd    = nobj/4;
  ipupil  = float(make_pupil(nobj,nobj/4,xc=cent,yc=cent,cobs=0.12));
  _p1     = long(ceil(cent-pupd/2.));
  _p2     = long(floor(cent+pupd/2.));
  
  prepzernike,nobj,pupd,cent,cent;
  modes = array(0.,pupd,pupd,nmodes);
  for (i=1;i<=100;i++) modes(,,i) = zernike(i+3)(_p1:_p2,_p1:_p2);

  npsf = 100;
  psfs = array(0.,nobj,nobj,npsf);
  ftos = array(complex,nobj,nobj,npsf);
  for (i=1;i<=npsf;i++) {
    coeffs    = random_n(nmodes)/8.;
    tmp       = mode2psf(nobj,ipupil(_p1:_p2,_p1:_p2),modes,coeffs,&phase,&amplipup,&amplifoc);
    tmp      /= sum(tmp);
    psfs(,,i) = tmp;
    ftos(,,i) = fft(tmp);
  }

  psfm = roll(psfs(,,avg));
  psfm = psfm/sum(psfm);
  ftom = fft(roll(psfm),1);
  dsp_psf = (abs(ftos-ftom)^2)(,,avg);
  clip,dsp_psf,max(dsp_psf)*1.e-6,;
  //dsp_psf(1) = 1;
  
  image=rfftconvol(object,psfs(,,npsf/2));
  psf1=roll(psfs(,,npsf/2));
  psf2=roll(psfs(,,1));

  varnoise = max(image/nlevel)^2;

  noise       = random_n([2,nobj,nobj])*sqrt(varnoise);
  image_noise = image+ noise;
  
  psf_noise   = (psf1*sum(image)+random_n([2,nobj,nobj])*sqrt(varnoise))/sum(image);

  write,"Estimating the object spectrum from the image";
  
  psdo=abs(fft(object,1))^2;
  psdo=clip(psdo,1.e-6,);
  //psdo=roll(psdo);
  psdo=psdo/max(psdo)*(sum(image_noise))^2;
  
  params=yoda_psdobj(image_noise,psf1,varn_est,psd_noise,psd_obj,disp=1,
                    verbose=0,nozero=0,threshold=1.e-6,noise_psf=1,object=object);

  ratio = psd_noise/psd_obj;

  psd_obj = psd_obj / max(psd_obj) * sum(image)^2;
  //computing wiener estimate
  myguess = fft((fft(image_noise)*fft(psf2)/(abs(fft(psf2))^2+ratio)),-1).re / numberof(image);
  myguess = roll(myguess);
  myguess /= sum(myguess);
  myguess *= sum(object);

  alpha = (max(ratio)>1. ? 1./max(ratio) : 1.);

  write,"Starting deconvolution with yoda ...";

  write,"Deconvolution with a gaussian prior";
  yoda,image=image_noise,psf=psf2,variance=varnoise,positivity=1,scale_psd = alpha ,psd=psdo,
    mobj=object,myopic = 1,scale_psf=1.,psd_psf = dsp_psf,disp=1,script=1,
    guess=_(clip(image_noise,1.e-6,),psf2),mpsf=psf1,nbiter=2000.;
  
  myopic = *yoda_data.object;
  myopic = myopic / avg(myopic) * avg(object);
  
  yoda,image=image_noise,psf=psf2,variance=varnoise,positivity=1,scale_psd= alpha ,psd=psdo,
    mobj=object,disp=1,script=1,guess=clip(image_noise,1.e-6,),nbiter=2000.;
  classical = *yoda_data.object;
  classical = classical / avg(classical) * avg(object);

  yoda,image=image_noise,psf=psf1,variance=varnoise,positivity=1,scale_psd= alpha ,psd=psdo,
    disp=1,script=1,guess=clip(image_noise,1.e-6,),nbiter=2000.;
  
  realpsf = *yoda_data.object;
  realpsf = realpsf / avg(realpsf) * avg(object);
  
  winkill,5;
  window,5,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = sqrt(abs(realpsf-object)^2);//realpsf-object;
  imrec(1:nobj,nobj+1:)  = sqrt(abs(classical-object)^2);//classical-object;
  imrec(nobj+1:,nobj+1:) = sqrt(abs(myopic-object)^2);//myopic-object;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Error (Myopic + Gaussian prior)";
  plt,"Class.",0.358,0.65,tosys=0,color="red",height=8;
  plt,"Myopic",0.395,0.65,tosys=0,color="red",height=8;
  plt,"Object",0.358,0.63,tosys=0,color="red",height=8;
  plt,"Real PSF",0.395,0.63,tosys=0,color="red",height=8;

  winkill,4;
  window,4,dpi=120;
  imrec = array(0.,2*nobj,2*nobj);
  imrec(1:nobj,1:nobj)   = object;
  imrec(nobj+1:,1:nobj)  = clip(image_noise,0.,);
  imrec(1:nobj,nobj+1:)  = classical;
  imrec(nobj+1:,nobj+1:) = myopic;
  pli,imrec;
  colorbar,min(imrec),max(imrec);
  pltitle,"Myopic deconvolution with Gaussian prior";
  plt,"Class.",0.358,0.65,tosys=0,color="red",height=8;
  plt,"Myopic",0.395,0.65,tosys=0,color="red",height=8;
  plt,"Object",0.358,0.63,tosys=0,color="red",height=8;
  plt,"Image",0.395,0.63,tosys=0,color="red",height=8;
  
  write,"wiener error     :",sqrt(sum(abs(myguess-object)^2)/numberof(object));
  write,"classical error  :",sqrt(sum(abs(classical-object)^2)/numberof(object));
  write,"myopic error     :",sqrt(sum(abs(myopic-object)^2)/numberof(object));
  write,"real psf error   :",sqrt(sum(abs(realpsf-object)^2)/numberof(object));
  return imrec;
}


/*
object=png_read("yoda2.png");
object=object(avg,,);
object=object(,::-1);
object(where(object <=63.75))=0;
object2=array(double,nobj,nobj);
object2(5:124,5:124)=object;
object=object2;
fits_write,"yoda.fits",object;
*/
