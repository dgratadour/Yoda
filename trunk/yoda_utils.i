
func mode2psf(size,mask,modes,coeff,&phase,&amplipup,&amplifoc,fto=)
/* DOCUMENT mode2psf(mask,modes,coeff,&phase,&amplipup,&amplifoc,fto=)
 *  
 * This routine returns the PSF (or the OTF) estimated from modes coefficients
 *
 */ 
{
  if (is_void(fto)) fto = 0;
  
  if (anyof(dimsof(mask)(2:3) != dimsof(modes)(2:3)))
    error,"size of mask and modes not compatible";

  pupd = dimsof(mask)(2);
  
  nbmodes = dimsof(modes)(4);
  
  if (nbmodes != numberof(coeff)) 
    error,"Sizes of modes and coeff not compatible";
  
  phase = coeff(+) * modes(,,+);
  
  psf=phase2psf(phase,mask,size,amplipup,amplifoc);
  
  if (fto) return fft(psf,1);
  else return psf;
}

func phase2psf(phase,mask,size,&amplipup,&amplifoc)
/* DOCUMENT phase2psf
 * psf=phase2psf(phase,mask,size,&amplipup,&amplifoc,fftws=)
 * This routine returns the PSF estimated from a phase map
 *
 * KEYWORDS :
 *
 * SEE ALSO phase2psfMode
 */ 
{
  pupd = dimsof(mask)(2);
  sphase = dimsof(phase);

  if (sphase(2) != pupd) 
    error,"Size of Phase and Mask not compatible";

  // Complex amplitude in the pupil plane
  amplipup = array(complex,[2,size,size]);
  (amplipup.re)(1:pupd,1:pupd) = cos(phase)*mask;
  (amplipup.im)(1:pupd,1:pupd) = sin(phase)*mask;

  // Complex amplitude in the focal plane
  amplifoc = fft(amplipup,1);
  psf      = abs(amplifoc)^2;

  // normalization
  psf      = psf/sum(psf);

  return psf;
}

func circavg(a,middle=)
/* DOCUMENT circavg(a,middle=)
 *  
 * returns the circular average of an array
 * if array is centered fourier-like (first pixel)
 * use middle flag
 *
 */ 
{
  dim=dimsof(a)(2);
  if (middle) return cart2pol(a)(,avg);
  else return cart2pol(roll(a))(,avg);
}
