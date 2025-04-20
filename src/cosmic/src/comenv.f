***
      SUBROUTINE COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &                  M02,M2,MC2,AJ2,JSPIN2,KW2,
     &                  ZPARS,ECC,SEP,JORB,COEL,star1,star2,vk,
     &                  kick_info,formation1,formation2,sigmahold,
     &                  bhspin1,bhspin2,binstate,mergertype,
     &                  jp,tphys,switchedCE,rad,tms,evolve_type,disrupt,
     &                  lumin,B_0,bacc,tacc,epoch,menv_bpp,renv_bpp,
     &                  bkick,deltam_1,deltam_2)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* Common Envelope Evolution.
*
*     Author : C. A. Tout
*     Date :   18th September 1996
*
*     Redone : J. R. Hurley
*     Date :   7th July 1998
*
*     Update : P. D. Kiel (for ECSN, fallback and bugs)
*     Date : cmc version mid 2010
*
*     Update : V. E. Delfavero (bugs in bpp reporting)
*     Date : 20th, April 2025
*
* Note on indexing binary components:
*
* M01, M1, MC1, AJ1, JSPIN1, KW1, formation1, bhspin1 are selected with star1
* M02, M2, MC2, AJ2, JSPIN2, KW2, formation2, bhspin1 are selected with star2
*
* deltam_1 and deltam_2 are calculated before being passed into comenv
*
* teff, radc, and ospin are calculated inside comenv (or a function)
*
* Other quantities are passed as vectors.
* These are lumin, menv_bpp, tms, rad, renv, bacc, tacc, epoch, and B_0.
* tms and rad are immediately separated into tms1/2_bpp and rad1/2_bpp
*
      INTEGER KW1,KW2,KW,KW1i,KW2i,snp
      INTEGER star1,star2
      INTEGER binstate,mergertype
*
      REAL*8 M01,M1,MC1,AJ1,JSPIN1,R1,L1,K21
      REAL*8 M02,M2,MC2,AJ2,JSPIN2,R2,L2,K22,MC22
      REAL*8 TSCLS1(20),TSCLS2(20),LUMS(10),GB(10),TM1,TM2,TN,ZPARS(20)
      REAL*8 EBINDI,EBINDF,EORBI,EORBF,ECIRC,SEPF,SEPL,MF,XX
      REAL*8 SEP_postCE, M_postCE
      REAL*8 CONST,DELY,DERI,DELMF,MC3,FAGE1,FAGE2
      REAL*8 ECC,SEP,JORB,TB,OORB,OSPIN1,OSPIN2,TWOPI
      REAL*8 RC1,RC2,Q1,Q2,RL1,RL2,LAMB1,LAMB2
      REAL*8 MENV,RENV,MENVD,RZAMS,vk
      REAL*8 Porbi,Porbf,Mcf,Menvf,qi,qf,G
      REAL*8 kick_info(2,18),fallback,M1i,M2i
      REAL*8 bkick(20)
      REAL*8 bhspin1,bhspin2
      REAL*8 deltam_1,deltam_2
      common /fall/fallback
      INTEGER formation1,formation2
      REAL*8 sigmahold
      REAL*8 AURSUN,K3
      PARAMETER (AURSUN = 214.95D0,K3 = 0.21D0)
      LOGICAL COEL,output
      REAL*8 CELAMF,RL,RZAMSF
      EXTERNAL CELAMF,RL,RZAMSF
      INTEGER jp
      REAL*8 tphys
      LOGICAL switchedCE,disrupt
      INTEGER kstar1_bpp,kstar2_bpp
      REAL*8 mass1_bpp,mass2_bpp
      REAL*8 rrl1_bpp,rrl2_bpp,evolve_type
      REAL*8 aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp
      REAL*8 massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp
      REAL*8 q1_bpp,q2_bpp,teff1,teff2
      REAL*8 KW1_TEMP, KW2_TEMP
      REAL*8 rad(2),tms(2),lumin(2),B_0(2),bacc(2),tacc(2),epoch(2)
      REAL*8 menv_bpp(2),renv_bpp(2)
*
* Initialize
*

      aj1_bpp = 0.d0
      aj2_bpp = 0.d0
      massc1_bpp = 0.d0
      massc2_bpp = 0.d0

      if(using_cmc.eq.0) then
         tms1_bpp = tms(1)
         tms2_bpp = tms(2)
         rad1_bpp = rad(1)
         rad2_bpp = rad(2)
      endif

*
* Common envelope evolution - entered only when KW1 = 2, 3, 4, 5, 6, 8 or 9.
*
* For simplicity energies are divided by -G.
*
      TWOPI = 2.D0*ACOS(-1.D0)
      COEL = .FALSE.
      snp = 0
      output = .false.
*
* Obtain the core masses and radii.
*
      KW = KW1
      CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
      CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &            R1,L1,KW1,MC1,RC1,MENV,RENV,K21,
     &            bhspin1,star1)
      OSPIN1 = JSPIN1/(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
      MENVD = MENV/(M1-MC1)
      RZAMS = RZAMSF(M01)
*
* Decide which CE prescription to use based on LAMBDA flag
* MJZ: NOTE - Nanjing lambda prescription DOES NOT WORK!
*
      LAMB1 = CELAMF(KW,M01,L1,R1,RZAMS,MENVD,LAMBDAF)
      KW = KW2
      CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
      CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &            R2,L2,KW2,MC2,RC2,MENV,RENV,K22,
     &            bhspin2,star2)
      OSPIN2 = JSPIN2/(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
*
* Calculate the binding energy of the giant envelope (multiplied by lambda).
*
      EBINDI = M1*(M1-MC1)/(LAMB1*R1)
*
* If the secondary star is also giant-like add its envelopes energy.
* Determine EORBI based on CEFLAG (CEFLAG=1 for de Kool prescription)
*
      IF(KW2.GE.2.AND.KW2.LE.9.AND.KW2.NE.7)THEN
         MENVD = MENV/(M2-MC2)
         RZAMS = RZAMSF(M02)
         LAMB2 = CELAMF(KW,M02,L2,R2,RZAMS,MENVD,LAMBDAF)
         EBINDI = EBINDI + M2*(M2-MC2)/(LAMB2*R2)
*
* Calculate the initial orbital energy
*
         IF(CEFLAG.EQ.0) EORBI = MC1*MC2/(2.D0*SEP)
         IF(CEFLAG.EQ.1) EORBI = M1*M2/(2.D0*SEP)
      ELSE
         IF(CEFLAG.EQ.0) EORBI = MC1*M2/(2.D0*SEP)
         IF(CEFLAG.EQ.1) EORBI = M1*M2/(2.D0*SEP)
      ENDIF
      if(output) write(*,*)'Init CE:',M01,M1,R1,M02,M2,R2,EBINDI,EORBI
*
* Allow for an eccentric orbit.
*
      ECIRC = EORBI/(1.D0 - ECC*ECC)
*
* Calculate the final orbital energy without coalescence.
*
      EORBF = EORBI + EBINDI/ALPHA1
*
* If the secondary is on the main sequence see if it fills its Roche lobe.
*
      IF(KW2.LE.1.OR.KW2.EQ.7)THEN
         SEPF = MC1*M2/(2.D0*EORBF)
         Q1 = MC1/M2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
*
* If cemergeflag is set, cause kstars without clear core-envelope
* structure to merge automatically if they enter a CE
*
         IF(cemergeflag.eq.1)THEN
            if(KW1.eq.0.or.KW1.eq.1.or.KW1.eq.2.or.KW1.eq.7.or.
     &         KW1.eq.8.or.KW1.eq.10.or.KW1.eq.11.or.KW1.eq.12)then
                  COEL = .TRUE.
                  SEPL = RC1/RL1
                  binstate = 1
                  CALL CONCATKSTARS(KW1, KW2, mergertype)
            endif
         ENDIF
*
* The helium core of a very massive star of type 4 may actually fill
* its Roche lobe in a wider orbit with a very low-mass secondary.
*
         IF(RC1/RL1.GE.R2/RL2)THEN
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
               binstate = 1
               CALL CONCATKSTARS(KW1, KW2, mergertype)
            ENDIF
         ELSE
            IF(R2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = R2/RL2
               binstate = 1
               CALL CONCATKSTARS(KW1, KW2, mergertype)
            ENDIF
         ENDIF
         IF(COEL)THEN
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1
            IF(KW2.EQ.7.AND.KW.EQ.4) MC3 = MC3 + M2
*
* Coalescence - calculate final binding energy.
*
            EORBF = MAX(MC1*M2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
            KW1_TEMP = KW
            KW2_TEMP = 15
         ELSE
*
* Primary becomes a black hole, neutron star, white dwarf or helium star.
*
            MF = M1
            M1 = MC1
            KW1i = KW1
            KW2i = KW2
            M1i = M1
*
* Choose which masses and separations to use in SN based on cekickflag
*
            IF(cekickflag.eq.0)then
               SEP_postCE=SEPF
               M_postCE=MF
            ELSEIF(cekickflag.eq.1)then
               SEP_postCE=SEP
               M_postCE=MF
            ELSEIF(cekickflag.eq.2)then
               SEP_postCE=SEPF
               M_postCE=MC1
            ENDIF

            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21,
     &                  bhspin1,star1)
            IF(KW1.GE.13)THEN
               formation1 = 1
               if(KW1.eq.13.and.ecsn.gt.0.d0)then
                  if(KW1i.le.6)then
                     if(M1i.le.zpars(5))then
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation1 = 2
                     endif
                  elseif(KW1i.ge.7.and.KW1i.le.9)then
                     if(M1i.gt.ecsn_mlow.and.M1i.le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation1 = 2
                     endif
                  elseif(formation1.eq.11)then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 5
                  elseif(KW1i.ge.10.or.KW1i.eq.12)then
* AIC formation, will never happen here but...
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 4
                  endif
               endif
*              Check if M1 and M2 were switched on pass to comenv
               if(switchedCE)then
                   mass1_bpp = M2
                   mass2_bpp = M_postCE
                   massc1_bpp = MC2
                   massc2_bpp = MC1

                   kstar1_bpp = KW2i
                   kstar2_bpp = KW1i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = R2/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = RC1/(RL(q2_bpp)*SEP_postCE)
                   aj1_bpp = AJ2
                   aj2_bpp = AJ1

                   evolve_type = 16.0d0

               else
                   mass1_bpp = M_postCE
                   mass2_bpp = M2
                   massc1_bpp = MC1
                   massc2_bpp = MC2

                   kstar1_bpp = KW1i
                   kstar2_bpp = KW2i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = RC1/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = R2/(RL(q2_bpp)*SEP_postCE)
                   evolve_type = 15.0d0
                   aj1_bpp = AJ1
                   aj2_bpp = AJ2
               endif
               TB = (SEP_postCE/AURSUN)*
     &               SQRT(SEP_postCE/(AURSUN*(M_postCE+M2)))
               if(using_cmc.eq.0)then
                   if(switchedCE)then
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0)) 
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M02,M01,lumin(1),lumin(2),teff1,teff2,
     &                      RC2,RC1,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN2,OSPIN1,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin2,bhspin1,
     &                      deltam_2,deltam_1,formation2,formation1,
     &                      binstate,mergertype,'bpp')
                   else
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M01,M02,lumin(1),lumin(2),teff1,teff2,
     &                      RC1,RC2,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN1,OSPIN2,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin1,bhspin2,
     &                      deltam_1,deltam_2,formation1,formation2,
     &                      binstate,mergertype,'bpp')
                   endif
               endif
               CALL kick(KW1,M_postCE,M1,M2,ECC,SEP_postCE,
     &                   JORB,vk,star1,R2,fallback,sigmahold,
     &                   kick_info,disrupt,bkick)
* Returning variable state to original naming convention
               MF = M_postCE
               SEPF = SEP_postCE
               snp = 1
               if(M2.lt.0.d0)then
                  if(KW2.ge.10) M1 = M1-M2
                  MC1 = M1
                  MC2 = 0.D0
                  bhspin2 = 0.D0
                  M2 = 0.D0
                  KW2 = 15
                  AJ1 = 0.D0
                  COEL = .true.
                  GOTO 30
               endif
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
         ENDIF
      ELSE
*
* Degenerate or giant secondary. Check if the least massive core fills its
* Roche lobe.
*
         SEPF = MC1*MC2/(2.D0*EORBF)
         Q1 = MC1/MC2
         Q2 = 1.D0/Q1
         RL1 = RL(Q1)
         RL2 = RL(Q2)
*
* If cemergeflag is set, cause kstars without clear core-envelope
* structure to merge automatically if they enter a CE
*
         IF(cemergeflag.eq.1)THEN
            if(KW1.eq.0.or.KW1.eq.1.or.KW1.eq.2.or.KW1.eq.7.or.
     &         KW1.eq.8.or.KW1.eq.10.or.KW1.eq.11.or.KW1.eq.12)then
                  COEL = .TRUE.
                  SEPL = RC1/RL1
                  binstate = 1
                  CALL CONCATKSTARS(KW1, KW2, mergertype)
            endif
         ENDIF
         IF(RC1/RL1.GE.RC2/RL2)THEN
            IF(RC1.GT.RL1*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC1/RL1
               binstate = 1
               CALL CONCATKSTARS(KW1, KW2, mergertype)
            ENDIF
         ELSE
            IF(RC2.GT.RL2*SEPF)THEN
               COEL = .TRUE.
               SEPL = RC2/RL2
               binstate = 1
               CALL CONCATKSTARS(KW1, KW2, mergertype)
            ENDIF
         ENDIF
*
         IF(COEL)THEN
*
* If the secondary was a neutron star or black hole the outcome
* is an unstable Thorne-Zytkow object that leaves only the core.
*
            SEPF = 0.D0
            IF(KW2.GE.13)THEN
               MC1 = MC2
               M1 = MC1
               MC2 = 0.D0
               M2 = 0.D0
               KW1 = KW2
               KW2 = 15
               bhspin1 = bhspin2
               bhspin2 = 0.d0
               AJ1 = 0.D0
*
* The envelope mass is not required in this case.
*
               GOTO 30
            ENDIF
*
            KW = KTYPE(KW1,KW2) - 100
            MC3 = MC1 + MC2
*
* Calculate the final envelope binding energy.
*
            EORBF = MAX(MC1*MC2/(2.D0*SEPL),EORBI)
            EBINDF = EBINDI - ALPHA1*(EORBF - EORBI)
            if(output) write(*,*)'In dg or giant 1:',M01,M1,R1,M02,M2,
     & R2,MC1,MC2,MC3,KW1,KW2,KW,EORBF,EBINDF
*
* Check if we have the merging of two degenerate cores and if so
* then see if the resulting core will survive or change form.
*
            IF(KW1.EQ.6.AND.(KW2.EQ.6.OR.KW2.GE.11))THEN
               CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
               if(output) write(*,*)'dg/giant 2:',KW,MC1,MC2,MC3,EBINDF
            ENDIF
            IF(KW1.LE.3.AND.M01.LE.ZPARS(2))THEN
               IF((KW2.GE.2.AND.KW2.LE.3.AND.M02.LE.ZPARS(2)).OR.
     &             KW2.EQ.10)THEN
                  CALL dgcore(KW1,KW2,KW,MC1,MC2,MC3,EBINDF)
                  if(output) write(*,*)'dg/giant 2:',KW,MC1,MC2,MC3,
     & EBINDF
                  IF(KW.GE.10)THEN
                     KW1 = KW
                     M1 = MC3
                     MC1 = MC3
                     IF(KW.LT.15) M01 = MC3
                     AJ1 = 0.D0
                     MC2 = 0.D0
                     M2 = 0.D0
                     KW2 = 15
                     binstate = 2
                     mergertype = -1
                     GOTO 30
                  ENDIF
               ENDIF
            ENDIF
*
         ELSE
*
* The cores do not coalesce - assign the correct masses and ages.
*
            MF = M1
            M1 = MC1
            KW1i = KW1
            KW2i = KW2
            M1i = M1
*
* Choose which masses and separations to use in SN based on cekickflag
*
            IF(cekickflag.eq.0)then
               SEP_postCE=SEPF
               M_postCE=MF
            ELSEIF(cekickflag.eq.1)then
               SEP_postCE=SEP
               M_postCE=MF
            ELSEIF(cekickflag.eq.2)then
               SEP_postCE=SEPF
               M_postCE=MC1
            ENDIF
*
* If cehestarflag is specified, determine the post-CE separation and/or
* mass for He-star/compact object system according to Tauris+2015
*
            IF(KW1.ge.7.and.KW1.le.9.and.KW2.ge.13)THEN
                if(cehestarflag.ne.0)then
* calculate G in Rsun/Msun/days, and initial Porb in days
                    G = 2957.45d+0
                    Porbi = (TWOPI/(G*(MF+M2))**(1d0/2d0)) *
     &                  SEP**(3d0/2d0)
                    Mcf = (1d0/(400d0*Porbi)+0.49d0)*MF -
     &                  ((0.016d0/Porbi) - 0.106d0)
                    if(Porbi.le.2d0.and.Porbi.gt.0.06d0)then
                        Menvf = 0.18d0*Porbi**(0.45d0) *
     &                      (LOG(Mcf**4d0) - 1.05d0)
                    elseif(Porbi.gt.2.d0)then
                        Menvf = Mcf*(LOG(Porbi**(-0.2d0))+1) +
     &                      LOG(Porbi**(0.5d0))-1.5d0
                    elseif(Porbi.le.0.06d0)then
*                       outside range of validity, use BSE masses
                        Menvf = 0d0
                        Mcf = MC1
                    endif
                    qi = MF/M2
* if cehestarflag is 1, use BSEs calculation of post-CE core mass
                    if(cehestarflag.eq.1)then
                        M_postCE = MC1
                        qf = MC1/M2
* elseif cehestarflag is 2, use Tauris fitting formula for post-CE mass
                    elseif(cehestarflag.eq.2)then
*                   first, check that the mass is more than core mass
                        if((Mcf+Menvf).lt.MC1)then
                           M_postCE = MC1
                        else
                           M_postCE = (Mcf+Menvf)
                        endif
                        qf = (Mcf+Menvf)/M2
                    endif
                    Porbf = (((qi+1)/(qf+1))**(2d0) * (qi/qf)**(3d0) *
     &                      EXP(3d0*(qf-qi))) * Porbi
                    SEP_postCE = (G*(M_postCE+M2) / (TWOPI**(2d0)) *
     &                      Porbf**(2d0))**(1d0/3d0)
                endif
            ENDIF

            CALL star(KW1,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &                  R1,L1,KW1,MC1,RC1,MENV,RENV,K21,
     &                  bhspin1,star1)
            IF(KW1.GE.13)THEN
               formation1 = 1
               if(KW1.eq.13.and.ecsn.gt.0.d0)then
                  if(KW1i.le.6)then
                     if(M1i.le.zpars(5))then
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation1 = 2
                     endif
                  elseif(KW1i.ge.7.and.KW1i.le.9)then
                     if(M1i.gt.ecsn_mlow.and.M1i.le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation1 = 2
                     endif
                  elseif(formation1.eq.11)then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 5
                  elseif(KW1i.ge.10.or.KW1i.eq.12)then
* AIC formation, will never happen here but...
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 4
                  endif
               endif
*              Check if M1 and M2 were switched on pass to comenv
               if(switchedCE)then
                   mass1_bpp = M2
                   mass2_bpp = M_postCE
                   massc1_bpp = MC2
                   massc2_bpp = MC1

                   kstar1_bpp = KW2i
                   kstar2_bpp = KW1i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = R2/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = RC1/(RL(q2_bpp)*SEP_postCE)
                   aj1_bpp = AJ2
                   aj2_bpp = AJ1
                   evolve_type = 16.d0

               else
                   mass1_bpp = M_postCE
                   mass2_bpp = M2
                   massc1_bpp = MC1
                   massc2_bpp = MC2

                   kstar1_bpp = KW1i
                   kstar2_bpp = KW2i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = RC1/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = R2/(RL(q2_bpp)*SEP_postCE)
                   aj1_bpp = AJ1
                   aj2_bpp = AJ2
                   evolve_type = 15.0d0
               endif
               TB = (SEP_postCE/AURSUN)*
     &               SQRT(SEP_postCE/(AURSUN*(M_postCE+M2)))
               if(using_cmc.eq.0)then
                   if(switchedCE)then
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M02,M01,lumin(1),lumin(2),teff1,teff2,
     &                      RC2,RC1,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN2,OSPIN1,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin2,bhspin1,
     &                      deltam_2,deltam_1,formation2,formation1,
     &                      binstate,mergertype,'bpp')
                   else
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M01,M02,lumin(1),lumin(2),teff1,teff2,
     &                      RC1,RC2,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN1,OSPIN2,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin1,bhspin2,
     &                      deltam_1,deltam_2,formation1,formation2,
     &                      binstate,mergertype,'bpp')
                   endif
               endif
* USSN: if ussn flag is set, have reduced kicks for stripped He stars (SN=8)
               if(KW1.eq.13.and.KW2.ge.13.and.ussn.eq.1)then
                  if(KW1i.ge.7.and.KW1i.le.9)then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                  formation1 = 3
                  endif
               endif
               CALL kick(KW1,M_postCE,M1,M2,ECC,SEP_postCE,
     &                   JORB,vk,star1,R2,fallback,sigmahold,
     &                   kick_info,disrupt,bkick)
* Returning variable state to original naming convention
               MF = M_postCE
               SEPF = SEP_postCE
               snp = 1
               if(M2.lt.0.d0)then
                  if(KW2.ge.10) M1 = M1-M2
                  MC1 = M1
                  MC2 = 0.D0
                  M2 = 0.D0
                  KW2 = 15
                  AJ1 = 0.D0
                  COEL = .true.
                  GOTO 30
               endif
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
            MF = M2
            KW = KW2
            M2 = MC2
            KW1i = KW1
            KW2i = KW2
            M2i = M2
*
* Choose which masses and separations to use in SN based on cekickflag
*
            IF(cekickflag.eq.0)then
               SEP_postCE=SEPF
               M_postCE=MF
            ELSEIF(cekickflag.eq.1)then
               SEP_postCE=SEP
               M_postCE=MF
            ELSEIF(cekickflag.eq.2)then
               SEP_postCE=SEPF
               M_postCE=MC2
            ENDIF

            CALL star(KW2,M02,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            CALL hrdiag(M02,AJ2,M2,TM2,TN,TSCLS2,LUMS,GB,ZPARS,
     &                  R2,L2,KW2,MC2,RC2,MENV,RENV,K22,
     &                  bhspin2,star2)
            IF(KW2.GE.13.AND.KW.LT.13)THEN
               formation2 = 1
               if(KW2.eq.13.and.ecsn.gt.0.d0)then
                  if(KW2i.le.6)then
                     if(M2i.le.zpars(5))then
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation2 = 2
                     endif
                  elseif(KW2i.ge.7.and.KW2i.le.9)then
                     if(M2i.gt.ecsn_mlow.and.M2i.le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                        formation2 = 2
                     endif
                  elseif(formation2.eq.11)then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation2 = 5
                  elseif(KW2i.ge.10.or.KW2i.eq.12)then
* AIC formation, will never happen here but...
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation2 = 4
                  endif
               endif
*              Check if M1 and M2 were switched on pass to comenv
               if(switchedCE)then
                   mass1_bpp = M_postCE
                   mass2_bpp = M1
                   massc1_bpp = MC2
                   massc2_bpp = MC1

                   kstar1_bpp = KW2i
                   kstar2_bpp = KW1i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = RC2/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = R1/(RL(q2_bpp)*SEP_postCE)
                   evolve_type = 15.0d0
                   aj1_bpp = AJ2
                   aj2_bpp = AJ1
*              M2, which here is mass(1), undergoes SN

               else
                   mass1_bpp = M1
                   mass2_bpp = M_postCE
                   massc1_bpp = MC1
                   massc2_bpp = MC2

                   kstar1_bpp = KW1i
                   kstar2_bpp = KW2i
                   q1_bpp = mass1_bpp/mass2_bpp
                   q2_bpp = 1.d0/q1_bpp
                   rrl1_bpp = R1/(RL(q1_bpp)*SEP_postCE)
                   rrl2_bpp = RC2/(RL(q2_bpp)*SEP_postCE)
                   evolve_type = 16.0d0
                   aj1_bpp = AJ1
                   aj2_bpp = AJ2
*              M2, which here is mass(2), undergoes SN
               endif
               TB = (SEP_postCE/AURSUN)*
     &               SQRT(SEP_postCE/(AURSUN*(M_postCE+M1)))
               if(using_cmc.eq.0)then
                   if(switchedCE)then
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M02,M01,lumin(1),lumin(2),teff1,teff2,
     &                      RC2,RC1,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN2,OSPIN1,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin2,bhspin1,
     &                      deltam_2,deltam_1,formation2,formation1,
     &                      binstate,mergertype,'bpp')
                   else
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,SEP_postCE,TB,ECC,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M01,M02,lumin(1),lumin(2),teff1,teff2,
     &                      RC1,RC2,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN1,OSPIN2,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin1,bhspin2,
     &                      deltam_1,deltam_2,formation1,formation2,
     &                      binstate,mergertype,'bpp')
                   endif
               endif
               CALL kick(KW2,M_postCE,M2,M1,ECC,SEP_postCE,
     &                   JORB,vk,star2,R1,fallback,sigmahold,
     &                   kick_info,disrupt,bkick)
* Returning variable state to original naming convention
               MF = M_postCE
               SEPF = SEP_postCE
               snp = 1
               if(M1.lt.0.d0)then
                  if(KW2.ge.10) M2 = M2-M1
                  MC2 = M2
                  MC1 = 0.D0
                  M1 = 0.D0
                  bhspin1 = 0.D0
                  KW1 = 15
                  AJ2 = 0.D0
                  COEL = .true.
                  GOTO 30
               endif
               IF(ECC.GT.1.D0) GOTO 30
            ENDIF
         ENDIF
      ENDIF
*
      IF(COEL)THEN
         MC22 = MC2
         IF(KW.EQ.4.OR.KW.EQ.7)THEN
* If making a helium burning star calculate the fractional age
* depending on the amount of helium that has burnt.
            IF(KW1.LE.3)THEN
               FAGE1 = 0.D0
            ELSEIF(KW1.GE.6)THEN
               FAGE1 = 1.D0
            ELSE
               FAGE1 = (AJ1 - TSCLS1(2))/(TSCLS1(13) - TSCLS1(2))
            ENDIF
            IF(KW2.LE.3.OR.KW2.EQ.10)THEN
               FAGE2 = 0.D0
            ELSEIF(KW2.EQ.7)THEN
               FAGE2 = AJ2/TM2
               MC22 = M2
            ELSEIF(KW2.GE.6)THEN
               FAGE2 = 1.D0
            ELSE
               FAGE2 = (AJ2 - TSCLS2(2))/(TSCLS2(13) - TSCLS2(2))
            ENDIF
         ENDIF
         if(output) write(*,*)'coel 1:',KW,MC1,MC2,MC3,MC22,FAGE1,FAGE2
      ENDIF
*
* Now calculate the final mass following coelescence.  This requires a
* Newton-Raphson iteration.
*
      IF(COEL)THEN
*
* Calculate the orbital spin just before coalescence.
*
         TB = (SEPL/AURSUN)*SQRT(SEPL/(AURSUN*(MC1+MC2)))
         OORB = TWOPI/TB
*
         XX = 1.D0 + ZPARS(7)
         IF(EBINDF.LE.0.D0)THEN
            MF = MC3
            GOTO 20
         ELSE
            CONST = ((M1+M2)**XX)*(M1-MC1+M2-MC22)*EBINDF/EBINDI
         ENDIF
*
* Initial Guess.
*
         MF = MAX(MC1 + MC22,(M1 + M2)*(EBINDF/EBINDI)**(1.D0/XX))
   10    DELY = (MF**XX)*(MF - MC1 - MC22) - CONST
*        IF(ABS(DELY/MF**(1.D0+XX)).LE.1.0D-02) GOTO 20
         IF(ABS(DELY/MF).LE.1.0D-03) GOTO 20
         DERI = MF**ZPARS(7)*((1.D0+XX)*MF - XX*(MC1 + MC22))
         DELMF = DELY/DERI
         MF = MF - DELMF
         GOTO 10
*
* Set the masses and separation.
*
   20    IF(MC22.EQ.0.D0) MF = MAX(MF,MC1+M2)
         M2 = 0.D0
         M1 = MF
         KW2 = 15
         bhspin2 = 0.d0
*
* Combine the core masses.
*
         if(output) write(*,*)'coel 2 1:',KW,KW1,KW2,M1,M2,MF,MC22,
     & TB,OORB
         IF(KW.EQ.2)THEN
            CALL star(KW,M1,M1,TM2,TN,TSCLS2,LUMS,GB,ZPARS)
            IF(GB(9).GE.MC1)THEN
               M01 = M1
               AJ1 = TM2 + (TSCLS2(1) - TM2)*(AJ1-TM1)/(TSCLS1(1) - TM1)
               CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            ENDIF
            if(output) write(*,*)'coel 2 2:',KW,KW1,KW2,M1,M01,MC22,
     & TB,OORB
         ELSEIF(KW.EQ.7)THEN
            M01 = M1
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            AJ1 = TM1*(FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
            if(output) write(*,*)'coel 2 3:',KW,KW1,KW2,M1,M01,MC22,
     & TB,OORB
         ELSEIF(KW.EQ.4.OR.MC2.GT.0.D0.OR.KW.NE.KW1)THEN
            IF(KW.EQ.4) AJ1 = (FAGE1*MC1 + FAGE2*MC22)/(MC1 + MC22)
            MC1 = MC1 + MC2
            MC2 = 0.D0
*
* Obtain a new age for the giant.
*
            CALL gntage(MC1,M1,KW,ZPARS,M01,AJ1)
            CALL star(KW,M01,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS)
            if(output) write(*,*)'coel 2 4:',KW,KW1,KW2,M1,M01,MC22,
     & TB,OORB
         ENDIF
         MF = M1
         KW1i = KW
         KW2i = KW2
         M1i = M1
         CALL hrdiag(M01,AJ1,M1,TM1,TN,TSCLS1,LUMS,GB,ZPARS,
     &               R1,L1,KW,MC1,RC1,MENV,RENV,K21,
     &               bhspin1,star1)
         if(output) write(*,*)'coel 2 5:',KW,M1,M01,R1,MENV,RENV
         IF(KW1i.LE.12.and.KW.GE.13)THEN
            formation1 = 1
            if(KW.eq.13.and.ecsn.gt.0.d0)then
               if(KW1i.le.6)then
                  if(M1i.le.zpars(5))then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 2
                  endif
               elseif(KW1i.ge.7.and.KW1i.le.9)then
                  if(M1i.gt.ecsn_mlow.and.M1i.le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation1 = 2
                  endif
               elseif(formation1.eq.11)then
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation1 = 5
               elseif(KW1i.ge.10.or.KW1i.eq.12)then
* AIC formation, will never happen here but...
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation1 = 4
               endif
            endif
*              Check if M1 and M2 were switched on pass to comenv
            if(switchedCE)then
                mass1_bpp = 0.d0
                mass2_bpp = MF
                massc1_bpp = 0.d0
                massc2_bpp = MC1

                kstar1_bpp = KW2i
                kstar2_bpp = KW1i
                rrl1_bpp = 0.d0
                rrl2_bpp = 0.d0
                evolve_type = 16.0d0
                aj1_bpp = AJ2
                aj2_bpp = AJ1
            else
                mass1_bpp = MF
                mass2_bpp = 0.d0
                massc1_bpp = MC1
                massc2_bpp = 0.d0

                kstar1_bpp = KW1i
                kstar2_bpp = KW2i
                rrl1_bpp = 0.d0
                rrl2_bpp = 0.d0
                evolve_type = 15.0d0
                aj1_bpp = AJ1
                aj2_bpp = AJ2
            endif
            TB = 0.d0
            if(using_cmc.eq.0)then
                   if(switchedCE)then
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,-1.d0,TB,0.d0,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M02,M01,lumin(1),lumin(2),teff1,teff2,
     &                      RC2,RC1,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN2,OSPIN1,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin2,bhspin1,
     &                      deltam_2,deltam_1,formation2,formation1,
     &                      binstate,mergertype,'bpp')
                   else
                       teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad1_bpp**2.d0))**(1.d0/4.d0))
                       teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad2_bpp**2.d0))**(1.d0/4.d0))
                       CALL writetab(jp,tphys,evolve_type,
     &                       mass1_bpp,mass2_bpp,kstar1_bpp,
     &                       kstar2_bpp,-1.d0,TB,0.d0,
     &                       rrl1_bpp,rrl2_bpp,
     &                       aj1_bpp,aj2_bpp,tms1_bpp,tms2_bpp,
     &                       massc1_bpp,massc2_bpp,rad1_bpp,rad2_bpp,
     &                      M01,M02,lumin(1),lumin(2),teff1,teff2,
     &                      RC1,RC2,menv_bpp(1),menv_bpp(2),renv_bpp(1),
     &                      renv_bpp(2),OSPIN1,OSPIN2,B_0(1),B_0(2),
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin1,bhspin2,
     &                      deltam_1,deltam_2,formation1,formation2,
     &                      binstate,mergertype,'bpp')
                   endif
            endif
            CALL kick(KW,MF,M1,0.d0,0.d0,-1.d0,0.d0,vk,star1,
     &                0.d0,fallback,sigmahold,kick_info,disrupt,bkick)
            if(output) write(*,*)'coel 2 6:',KW,M1,M01,R1,MENV,RENV
         ENDIF
         JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         KW1 = KW
         ECC = 0.D0
         if(output) write(*,*)'coel 2 7:',KW1,M1,M01,R1,MENV,RENV
      ELSE
*
* Check if any eccentricity remains in the orbit by first using
* energy to circularise the orbit before removing angular momentum.
* (note this should not be done in case of CE SN ... fixed PDK).
*
         IF(snp.EQ.0)THEN
            IF(EORBF.LT.ECIRC)THEN
               ECC = SQRT(1.D0 - EORBF/ECIRC)
            ELSE
               ECC = 0.D0
            ENDIF
         ENDIF
*
* Set both cores in co-rotation with the orbit on exit of CE,
*
         TB = (SEPF/AURSUN)*SQRT(SEPF/(AURSUN*(M1+M2)))
         OORB = TWOPI/TB
         JORB = M1*M2/(M1+M2)*SQRT(1.D0-ECC*ECC)*SEPF*SEPF*OORB
*        JSPIN1 = OORB*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
*        JSPIN2 = OORB*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
*
* or, leave the spins of the cores as they were on entry.
* Tides will deal with any synchronization later.
*
         JSPIN1 = OSPIN1*(K21*R1*R1*(M1-MC1)+K3*RC1*RC1*MC1)
         JSPIN2 = OSPIN2*(K22*R2*R2*(M2-MC2)+K3*RC2*RC2*MC2)
      ENDIF
   30 SEP = SEPF
      if(output) write(*,*)'end of CE1:',KW1,M1,M01,R1,MENV,RENV
      if(output) write(*,*)'end of CE1:',KW2,M2,M02,R2,MENV,RENV
      sigma = sigmahold
      RETURN
      END
***
