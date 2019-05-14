*
* const_bse.h
*
      INTEGER idum1
      COMMON /RAND1/ idum1
      INTEGER idum2,iy,ir(32)
      COMMON /RAND2/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER tflag,ifflag,nsflag,wdflag,bhflag,windflag,ppsn
      COMMON /FLAGS/ tflag,ifflag,nsflag,wdflag,bhflag,windflag,ppsn
      INTEGER ceflag,cekickflag,cemergeflag,cehestarflag,ussn
      COMMON /CEFLAGS/ ceflag,cekickflag,cemergeflag,cehestarflag,ussn
*
      REAL*8 neta,bwind,hewind,mxns,beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /WINDVARS/ neta,bwind,hewind,mxns,beta,xi,acc2,epsnov
      COMMON /WINDVARS/ eddfac,gamma
      REAL*8 alpha1,lambda
      COMMON /CEVARS/ alpha1,lambda
      REAL*8 bconst,CK
      COMMON /MAGVARS/ bconst,CK
      REAL*8 sigma,bhsigmafrac,polar_kick_angle,mu_SN1,omega_SN1
      COMMON /SNVARS/ sigma,bhsigmafrac,polar_kick_angle
      COMMON /SNVARS/ mu_SN1,omega_SN1
*
      INTEGER*8 id1_pass,id2_pass
      REAL*8 merger
      COMMON /CMCPASS/ merger, id1_pass,id2_pass
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL*8 scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL*8 bcm(50000,42),bpp(1000,15)
      COMMON /BINARY/ bcm,bpp
*
