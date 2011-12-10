#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
#include "/work/home/eric/work/deconvolution/yan-0.1.0/deconv_sn1987a_init.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1)
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1)
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1)

info,lkl_h
dbexit
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1)
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/optimpack.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0)
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0)

m
n
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
dbexit
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0)
#include "/home/eric/work/yeti-5.3.3/optimpack/yorick/lbfgsb.i"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0,beta=0)
pli,x;
palette,"stern.gp"
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=0,beta=0,mu=0.1)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=0.1)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=0.1)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=10.)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=1,mu=10.)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=1,mu=10.)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=1,mu=10.)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=1,mu=1e2)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=1,mu=1e4)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e4)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e6)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e9)
fma;pli,x;
x=[]
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e9)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e9)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e9)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e9,frtol=1e-15)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e8,frtol=1e-15)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e7,frtol=1e-15)
fma;pli,x;
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e6,frtol=1e-15)
x = op_deconv(img,psf,x,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e6,frtol=1e-15)
fma;pli,x;
fma;pli,log(x);
x = op_deconv(img,psf,,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e6,frtol=1e-15)
fma;pli,x;
fma;pli,log(x);
x = op_deconv(img,psf,,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e6,frtol=1e-15)
#include "/work/home/eric/work/yeti-5.3.3/optimpack/yorick/op_deconv.i"
x = op_deconv(img,psf,,entropy=1,maxiter=10,verb=1,method=1,beta=0,mu=1e6,frtol=1e-15)
x = op_deconv(img,psf,,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e6,frtol=1e-15)
x = op_deconv(img,psf,,entropy=1,maxiter=10,verb=1,method=2,beta=0,mu=1e7,frtol=1e-15)
fma;pli,x;
