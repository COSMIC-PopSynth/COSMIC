***
      PROGRAM bse
***
*
* Evolves a binary by calling evolv2.f
* (see header of subroutine for algorithm description).
*
* This is a debugging driver: it evolves a single binary in the same way
* that cosmic.evolve does so that it can be stepped through with gdb.
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
      integer i,j,k,kstar(2)
      integer bpp_index_out, bcm_index_out
*
      real*8 mass0(2),mass(2),z,zpars(20),kick_info(2,19)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc
      real*8 B_0(2),bacc(2),tacc(2),bhspin(2)
*
************************************************************************
* Input:
*
* binary.in is written by debug/create_binary_in.py from a row of an
* initC table and contains three namelists, which are matched by name
* (so the order within each one does not matter):
*
*   &binary_in          the initial state of the binary, i.e. everything
*                       that cosmic.evolve passes as arguments to evolv2
*   &bse_in             the BSE settings, i.e. everything in
*                       GROUPED_SETTINGS in consts.py (randomseed=idum1)
*   &stellar_engine_in  which stellar engine to use and its settings
*
* Anything missing from binary.in keeps the default set below. Anything
* in binary.in that is not in these namelists stops the read with an
* error, so if you add a setting to COSMIC, add it to the namelist here.
*
* mass is in solar units.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state.
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* If using_metisse = 1 then the METISSE tracks are read from
* metisse_tracks.bin, which create_binary_in.py writes alongside binary.in
* (see test_metisse_tracks.f90).
*
      NAMELIST /binary_in/ kstar,mass0,mass,tb,ecc,z,tphysf,dtp,
     &         rad,lum,massc,radc,menv,renv,ospin,B_0,bacc,tacc,
     &         epoch,tms,bhspin,tphys
      NAMELIST /bse_in/ neta,bwind,hewind,beta,xi,acc2,epsnov,
     &         eddfac,gamma,LBV_flag,
     &         alpha1,lambdaf,qcrit_array,
     &         ceflag,cekickflag,cemergeflag,cehestarflag,ussn,
     &         tflag,ifflag,wdflag,rtmsflag,bhflag,remnantflag,
     &         maltsev_mode,grflag,bhms_coll_flag,wd_mass_lim,aic,
     &         bdecayfac,windflag,qcflag,eddlimflag,bhspinflag,
     &         rejuvflag,htpmb,ST_cr,ST_tide,
     &         pisn,ppi_co_shift,ppi_extra_ml,maltsev_fallback,
     &         maltsev_pf_prob,fryer_mass_limit,mxns,fryer_fmix,
     &         fryer_mcrit_nsbh,ecsn,ecsn_mlow,sigma,sigmadiv,
     &         bhsigmafrac,polar_kick_angle,natal_kick_array,
     &         bhspinmag,rembar_massloss,kickflag,mm_mu_ns,mm_mu_bh,
     &         pts1,pts2,pts3,
     &         don_lim,acc_lim,smt_periastron_check,
     &         bconst,CK,
     &         fprimc_array,
     &         rejuv_fac,
     &         zsun,
     &         idum1
      NAMELIST /stellar_engine_in/ using_sse,using_metisse,
     &         path_to_tracks,path_to_he_tracks,z_match_limit,
     &         METISSE_verbose
*
* Defaults for the initial state: start both stars on the ZAMS at t=0.
* If ospin is zero (actually < 0.001) at time zero then evolv2 will set
* an appropriate ZAMS spin.
*
      kstar = 1
      mass0 = 0.d0
      mass = -1.d0
      tb = 0.d0
      ecc = 0.d0
      z = 0.02d0
      tphysf = 13700.d0
      dtp = 0.d0
      rad = 0.d0
      lum = 0.d0
      massc = 0.d0
      radc = 0.d0
      menv = 0.d0
      renv = 0.d0
      ospin = 0.d0
      B_0 = 0.d0
      bacc = 0.d0
      tacc = 0.d0
      epoch = 0.d0
      tms = 0.d0
      bhspin = 0.d0
      tphys = 0.d0
*
      using_sse = 1
      using_metisse = 0
      path_to_tracks = ''
      path_to_he_tracks = ''
      z_match_limit = 0.d0
      METISSE_verbose = .false.
*
      OPEN(22,file='binary.in', status='old')
      READ(22,NML=binary_in)
      READ(22,NML=bse_in)
      READ(22,NML=stellar_engine_in)
      CLOSE(22)
*
      if(mass(1).lt.0.d0) mass(1) = mass0(1)
      if(mass(2).lt.0.d0) mass(2) = mass0(2)
*
      WRITE(*,*)'mass0, tphysf, tb, kstar, z, ecc:'
      WRITE(*,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      if(using_metisse.eq.1)then
         WRITE(*,*)'Using METISSE with tracks from:'
         WRITE(*,*)'  ',TRIM(path_to_tracks)
         WRITE(*,*)'  ',TRIM(path_to_he_tracks)
      endif
*
* Mirror what cosmic.evolve does before calling evolv2.
* If a star is already a compact object then one SN has already occurred.
*
      kick_info = 0.d0
      if((kstar(1).ge.10.and.kstar(1).le.14).or.
     &   (kstar(2).ge.10.and.kstar(2).le.14)) kick_info(1,1) = 1.d0
      using_cmc = 0
*
* Save every column (in the order of ALL_COLUMNS in consts.py) in bpp/bcm.
*
      n_col_bpp = 52
      n_col_bcm = 52
      do i = 1,52
         col_inds_bpp(i) = i
         col_inds_bcm(i) = i
      enddo
*
* Pass the METISSE tracks to the fortran in the same way python does.
* evolv2 itself then calls zcnsts (and instar) so we don't do that here.
* zpars must be zero so that METISSE processes the tracks from scratch, as
* it does for the first binary evolved in a python session.
*
      if(using_metisse.eq.1)then
         call load_metisse_tracks('metisse_tracks.bin')
      endif
      zpars = 0.d0
      WRITE(*,*)
*
* Evolve the binary.
*
      CALL evolv2(kstar,mass,tb,ecc,z,tphysf,dtp,mass0,rad,lum,
     &            massc,radc,menv,renv,ospin,B_0,bacc,tacc,
     &            epoch,tms,bhspin,tphys,zpars,kick_info,
     &            bpp_index_out,bcm_index_out)
*
************************************************************************
* Output:
*
* bpp.dat and bcm.dat contain the bpp and bcm tables with one row per line
* and the columns in the order of ALL_COLUMNS in consts.py, e.g. read with
*   pd.read_csv('bpp.dat', sep=r'\s+', names=consts.ALL_COLUMNS)
*
      if(bpp_index_out.lt.0)then
         WRITE(*,*)'evolv2 failed (METISSE error in zcnsts?)'
         STOP 1
      endif
*
      OPEN(23,file='bpp.dat', status='unknown')
      do j = 1,bpp_index_out
         WRITE(23,99)(bpp(j,k),k=1,52)
      enddo
      CLOSE(23)
*
      OPEN(23,file='bcm.dat', status='unknown')
      do j = 1,bcm_index_out
         WRITE(23,99)(bcm(j,k),k=1,52)
      enddo
      CLOSE(23)
 99   FORMAT(52es26.16e3)
*
      WRITE(*,*)'Wrote',bpp_index_out,' rows to bpp.dat and',
     &          bcm_index_out,' rows to bcm.dat'
*
************************************************************************
*
      STOP
      END
***
