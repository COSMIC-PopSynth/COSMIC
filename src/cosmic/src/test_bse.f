***
      PROGRAM bse
***
*
* Evolves a binary by calling evolv2.f 
* (see header of subroutine for algorithm description). 
*
* Required input is described below. 
***
* See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
* Updated reference is:
*           Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
* (please use this one).
***
* For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
* or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
* The binary evolution algorithm is described in Chapter 3 of the thesis.
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
      implicit none
*
      INCLUDE 'const_bse.h'
*
      integer i,kw,kw2,kstar(2),j,k,time,idum
      integer bpp_index_out, bcm_index_out
      real*8 kick_info_out(2,17)
*
      real*8 mass0(2),mass(2),z,zpars(20),vs(3),kick_info(2,17)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc,yearsc
      real*8 B_0(2),bacc(2),tacc(2),bhspin(2),bkick(20)
      PARAMETER(yearsc=3.1557d+07)
      CHARACTER*8 label(16)
*
************************************************************************
* Input:
*
* mass is in solar units.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag = 0 activates spin-energy correction in common-envelope (0). 
* ceflag = 1 activates de Kool common-envelope model (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* rtmsflag > 0 uses rtms from simulation data (1=[Boost], 2=[Bpass]) (0)
* remnantflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, remnantflag=0; 3.0, remnantflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
* If you enter a negative kstar then parameters for an evolved star are
* required in the order of:
* current age, initial mass and spin rate, 
* otherwise the star will start on the ZAMS.
*
      OPEN(22,file='binary.in', status='old')
      READ(22,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      READ(22,*)neta,bwind,hewind,alpha1,lambdaf,windflag,rtmsflag
      READ(22,*)ceflag,tflag,ifflag,wdflag,bhflag,remnantflag,mxns,idum
      READ(22,*)pts1,pts2,pts3
      READ(22,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
      if(kstar(1).lt.0.or.kstar(2).lt.0)then
         READ(22,*)tphys
         READ(22,*)aj,mass(1),ospin(1)
         epoch(1) = tphys - aj
         kstar(1) = ABS(kstar(1))
         READ(22,*)aj,mass(2),ospin(2)
         epoch(2) = tphys - aj
         kstar(2) = ABS(kstar(2))
      else

      WRITE(*,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      WRITE(*,*)neta,bwind,hewind,alpha1,lambdaf,windflag,rtmsflag
      WRITE(*,*)ceflag,tflag,ifflag,wdflag,bhflag,remnantflag,mxns,idum
      WRITE(*,*)pts1,pts2,pts3
      WRITE(*,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
*
* Initialize the parameters.
* Set the initial spin of the stars. If ospin is zero (actually < 0.001)
* at time zero then evolv2 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
* If ospin is negative then the stars will be in co-rotation with the orbit.
*
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
         epoch(2) = 0.d0
         ospin(2) = 0.d0

         B_0(1) = 0.d0
         B_0(2) = 0.d0
         bacc(1) = 0.d0
         bacc(2) = 0.d0
         tacc(1) = 0.d0
         tacc(2) = 0.d0
         bhspin(1) = 0.d0
         bhspin(2) = 0.d0
         tms(1) = 0.d0
         tms(2) = 0.d0
*
* NOTE: The following parameters are set to current default values (v3.4.10).
* They should really be included in binary.in and copied from an initC
* file. This is a temporary measure until the input file is updated.
*
         pisn = 45.d0
         cekickflag = 2
         cehestarflag = 0
         grflag = 1
         bhms_coll_flag = 0
         wd_mass_lim = 1
         ecsn = 2.25
         ecsn_mlow = 1.6
         aic = 1
         ussn = 0
         sigmadiv = -20
         bhsigmafrac = 1.0
         polar_kick_angle = 90
         do i = 1,4
            do j = 1,2
               natal_kick_array(j, i) = -100.d0
            enddo
        enddo
        natal_kick_array(1,5) = 0.0
        natal_kick_array(2,5) = 0.0
        do i = 1,8
            qcrit_array(i) = 0.0
        enddo
        don_lim = -1
        acc_lim = -1
        bdecayfac = 1
        bconst = 3000
        ck = 1000
        qcflag = 5
        eddlimflag = 0
        do i = 1,16
            fprimc_array(i) = 2.0/21.0
        enddo
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        st_cr = 1
        st_tide = 1
        rembar_massloss = 0.5
        zsun = 0.014
        kickflag = 0
        using_cmc = 0


         do i=1,20
            bkick(i) = 0.d0
         enddo 
         

      endif
      do i=1,17
        do j=1,2
          kick_info(j,i) = 0.d0
        enddo
      enddo
* If you would like to enter the seperation as input in place of the binary
* orbital period uncomment these lines (depending upon which units you wish
* to use).
*      tb = sqrt((tb/aursun)**3/(mass0(1)+mass0(2))) !if input was separation [Rsun] use this line.
*      tb = sqrt((tb)**3/(mass0(1)+mass0(2))) !if input was separation [AU] use this line.
*      tb = tb*yeardy
      if(idum.gt.0) idum = -idum
      CLOSE(22)
      WRITE(*,*)
*
* Note that this routine can be used to evolve a single star if you 
* simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
* well as some dummy value for ecc). 
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
*
* Set the collision matrix.
*
      CALL instar
*
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'
      label(15) = 'SN PRIM'
      label(16) = 'SN SEC'
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0
*
* Evolve the binary.
* 
      CALL evolv2(kstar,mass,tb,ecc,z,tphysf,dtp,mass0,rad,lum,
     &            massc,radc,menv,renv,ospin,B_0,bacc,tacc,
     &            epoch,tms,bhspin,tphys,zpars,bkick,kick_info,
     &            bcm_index_out,bpp_index_out,kick_info_out)
*
************************************************************************
* Output:
* First check that bcm is not empty.
*
      if(bcm(1,1).lt.0.0) goto 50
*
* The bcm array stores the stellar and orbital parameters at the 
* specified output times. The parameters are (in order of storage):
*
*    Time, 
*    [stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, core radius, mass of any convective 
*    envelope, radius of the envelope, epoch, spin, mass loss rate and 
*    ratio of radius to roche lobe radius (repeated for secondary)],
*    period, separation, eccentricity.
*
      OPEN(23,file='binary.dat', status='unknown')
      j = 0
 30   j = j + 1
      if(bcm(j,1).lt.0.0)then
         bcm(j-1,1) = bcm(j,1)
         j = j - 1
      endif
      kw = INT(bcm(j,2))
      kw2 = INT(bcm(j,16))
      WRITE(23,99)bcm(j,1),kw,kw2,bcm(j,4),bcm(j,18),
     &            bcm(j,8),bcm(j,22), 
     &            bcm(j,6),bcm(j,20),bcm(j,15),bcm(j,29),
     &            bcm(j,5),bcm(j,19),bcm(j,13),bcm(j,27),
     &            bcm(j,14),bcm(j,28),
     &            bcm(j,31),bcm(j,32)
      if(bcm(j,1).ge.0.0) goto 30
      CLOSE(23)
 99   FORMAT(g30.18,2i3,10g30.18,5e30.18,g30.18)
 999  FORMAT(g30.18,2g30.18,1p,2e30.18)
*
* The bpp array acts as a log, storing parameters at each change
* of evolution stage.
*
 50   j = 0
      WRITE(*,*)'     TIME      M1       M2   K1 K2        SEP    ECC',  
     &          '  R1/ROL1 R2/ROL2  TYPE'
 52   j = j + 1
      if(bpp(j,1).lt.0.0) goto 60
      kstar(1) = INT(bpp(j,4))
      kstar(2) = INT(bpp(j,5))
      kw = INT(bpp(j,11))
      WRITE(*,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw)
      goto 52
 60   continue
 100  FORMAT(g30.18,2g30.18,2i3,g30.18,g30.18,2g30.18,2x,a8)
      WRITE(*,*)
*
************************************************************************
*
      WRITE(*,*) bcm(2,31), bpp(j, 6)

      STOP
      END
***
