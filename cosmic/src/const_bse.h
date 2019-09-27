*
* const_bse.h
*
      INTEGER idum1
      COMMON /RAND1/ idum1
      INTEGER idum2,iy,ir(32)
      COMMON /RAND2/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER tflag,ifflag,nsflag,wdflag,bhflag,windflag,qcflag
      INTEGER eddlimflag,bhspinflag,aic
      COMMON /FLAGS/ tflag,ifflag,nsflag,wdflag,bhflag,windflag,qcflag
      COMMON /FLAGS/ eddlimflag,bhspinflag,aic
      INTEGER ceflag,cekickflag,cemergeflag,cehestarflag,ussn
      COMMON /CEFLAGS/ ceflag,cekickflag,cemergeflag,cehestarflag,ussn
      INTEGER pisn_track(2)
      COMMON /TRACKERS/ pisn_track
*
      REAL*8 neta,bwind,hewind,mxns,beta,xi,acc2,epsnov
      REAL*8 eddfac,gamma
      COMMON /WINDVARS/ neta,bwind,hewind,mxns,beta,xi,acc2,epsnov
      COMMON /WINDVARS/ eddfac,gamma
      REAL*8 alpha1,lambdaf
      REAL*8 qcrit_array(16)
      COMMON /CEVARS/ qcrit_array,alpha1,lambdaf
      REAL*8 bconst,CK
      COMMON /MAGVARS/ bconst,CK
      REAL*8 sigma,sigmadiv,bhsigmafrac,pisn
      REAL*8 polar_kick_angle,mu_SN1,omega_SN1
      REAL*8 ecsn,ecsn_mlow,bhspinmag
      REAL*8 natal_kick_array(6)
      COMMON /SNVARS/ natal_kick_array,sigma,sigmadiv,bhsigmafrac
      COMMON /SNVARS/ polar_kick_angle,mu_SN1,omega_SN1
      COMMON /SNVARS/ pisn,ecsn,ecsn_mlow,bhspinmag
      REAL*8 fprimc_array(16)
      COMMON /TIDALVARS/ fprimc_array
      REAL*8 rejuv_fac
      COMMON /MIXVARS/ rejuv_fac
*
      INTEGER*8 id1_pass,id2_pass,using_cmc
      REAL*8 merger
      COMMON /CMCPASS/ merger,id1_pass,id2_pass,using_cmc
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL*8 scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL*8 bcm(50000,42),bpp(1000,23)
      COMMON /BINARY/ bcm,bpp
*
