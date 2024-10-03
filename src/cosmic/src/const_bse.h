*
* const_bse.h
*
      INTEGER idum1
      COMMON /RAND1/ idum1
      INTEGER idum2,iy,ir(32)
      COMMON /RAND2/ idum2,iy,ir
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER tflag,ifflag,remnantflag,wdflag,bhflag,windflag,qcflag
      INTEGER eddlimflag,bhspinflag,aic,rejuvflag,rtmsflag
      INTEGER htpmb,ST_cr,ST_tide,bdecayfac,grflag,bhms_coll_flag
      INTEGER wd_mass_lim
      COMMON /FLAGS/ tflag,ifflag,remnantflag,wdflag,bhflag,windflag,
     &               qcflag,eddlimflag,bhspinflag,aic,rejuvflag,
     &               htpmb,ST_cr,ST_tide,bdecayfac,grflag,
     &               bhms_coll_flag,wd_mass_lim,rtmsflag
      REAL*8 don_lim,acc_lim,Mbh_initial
      COMMON /MTVARS/ don_lim,acc_lim,Mbh_initial
      INTEGER ceflag,cekickflag,cemergeflag,cehestarflag,ussn
      COMMON /CEFLAGS/ ceflag,cekickflag,cemergeflag,cehestarflag,ussn
      INTEGER pisn_track(2)
      COMMON /TRACKERS/ pisn_track
*
      REAL*8 zsun
      COMMON /METVARS/ zsun
      REAL*8 neta,bwind,hewind,beta,xi,acc2,epsnov
      REAL*8 eddfac,gamma
      COMMON /WINDVARS/ neta,bwind,hewind,beta,xi,acc2,epsnov,
     &                  eddfac,gamma
      REAL*8 alpha1,lambdaf
      REAL*8 qcrit_array(16)
      COMMON /CEVARS/ qcrit_array,alpha1,lambdaf
      REAL*8 bconst,CK
      COMMON /MAGVARS/ bconst,CK
      INTEGER kickflag
      REAL*8 sigma,sigmadiv,bhsigmafrac,pisn,mxns
      REAL*8 polar_kick_angle
      REAL*8 ecsn,ecsn_mlow,bhspinmag,rembar_massloss
      REAL*8 natal_kick_array(2,5)
      COMMON /SNVARS/ natal_kick_array,sigma,sigmadiv,bhsigmafrac,
     &            polar_kick_angle,pisn,ecsn,ecsn_mlow,
     &            bhspinmag,mxns,rembar_massloss,kickflag
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
      REAL*8 bcm(50000,49),bpp(1000,49)
      COMMON /BINARY/ bcm,bpp
      INTEGER n_col_bpp, n_col_bcm
      INTEGER col_inds_bpp(49), col_inds_bcm(49)
      COMMON /COL/ n_col_bpp,col_inds_bpp,n_col_bcm,col_inds_bcm
*
