***
      SUBROUTINE evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
     &                  epoch,tm,tphys,tphysf,dtp,z,zpars,kick_info)
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c     Plots the HRD and variables as a function of time.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
      INCLUDE 'const_bse.h'
*
      integer kw,it,ip,jp,j,kwold,fb,k2str,pulsar
      integer nv,intpol,k
      parameter(nv=50000)
*
      real*8 mass,z,aj
      real*8 epoch,tphys,tphys2,tmold,tbgold
      real*8 mt,tm,tn,tphysf,dtp,tsave
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,teff,rc,menv,renv,kick_info(2,17)
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
      real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
      real*8 tol,tiny
      parameter(tol=1.0d-10,tiny=1.0d-14)
      real*8 ajhold,rm0,eps,alpha2
      parameter(eps=1.0d-06,alpha2=0.09d0)
      real*8 mlwind,vrotf
      external mlwind,vrotf
      logical iplot,isave,disrupt
      common /fall/fallback
      REAL*8 kw3,wsun,wx
      PARAMETER(kw3=619.2d0,wsun=9.46d+07,wx=9.46d+08)
      REAL*8 yeardy,yearsc,aursun
      PARAMETER(yeardy=365.24d0,aursun=214.95d0,yearsc=3.1557d+07)
      REAL*8 twopi,u1,u2,mew,tacc
      REAL*8 fallback,vk,B_0,B,ospbru,b01_bcm,b_mdot,b_mdot_lim
      REAL*8 bacc,Bbot,convradcomp,deltam,dspint,eqspin
      REAL*8 dtj,evolve_type
      REAL*8 omdot,sigmahold,sn,bhspin,s
      REAL*8 djtx,djspint,jspbru,Kconst
*
      REAL ran3
      EXTERNAL ran3
*
      dtm = 0.d0
      r = 0.d0
      lum = 0.d0
      mc = 0.d0
      mc1 = 0.d0
      rc = 0.d0
      rl = 0.d0
      if(ospin.le.0.d0)then
         ospin = 1.0d-10
         jspin = 1.0d-10
      endif
      k2 = 0.15d0

* PDK additions
      fb = 1 !turns kick velocity limit owning to fallback on (1) or off (0).
      vk = 0.d0 !gets passed to and from kick(), is kick mag. can be used to set initial pulsar particulars.
      B_0 = 0.0
      k2str = 0
      if(kw.ne.13)then
            bacc = 0.d0
            tacc = 0.d0
      endif
      aic = 1
      b_mdot_lim = -1.0e-11
      Bbot = 5e+7
      Kconst = 2.5d-49

      twopi = 2.d0*ACOS(-1.d0)

*
* Setup variables which control the output (if it is required).
*
      ip = 1
      jp = 1
      tsave = tphys
      isave = .true.
      iplot = .false.
      if(dtp.le.0.d0)then
         iplot = .true.
         isave = .false.
         tsave = tphysf
      elseif(dtp.gt.tphysf)then
         isave = .false.
         tsave = tphysf
      endif
* 
      do 10 , j = 1,nv
*
         if(neta.gt.tiny.and.j.gt.1)then
*
* Calculate mass loss from the previous timestep.
*
            dt = 1.0d+06*dtm
            dms = mlwind(kw,lum,r,mt,mc,rl,z,tphys)*dt
            if(kw.lt.10)then
               dml = mt - mc
               if(dml.lt.dms)then
                  dtm = (dml/dms)*dtm
                  dms = dml
               endif
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Evaluate convective/radiative limits for a variety of stars as based
* on the work of Belczynski et al. (2008). As option.
* Note only certain stars of type k = 0, 1, 2, 4, 5, 6, 9 **double check
* 0, 1 have range of 0.35-Mms,conv. Mms,conv is function of metalicity.
* 2 & 4 have a temperature dependence. Convective if Teff < 10**3.73
* 3, 5 & 6 are giants and have no radiative envelopes.
* 9's envelope is convective when M < Mhe,conv, Mhe,conv = 3.0Msun.
*
         if(kw.le.1)then
*           Main sequence mass limit varying with metallicity
            convradcomp = 1.25d0
            if(z.gt.0.001d0.and.z.lt.zsun)then
               convradcomp = 0.747d0 + 55.73d0*z - 1532*z*z
            elseif(z.le.0.001d0)then
               convradcomp = 0.8d0
            endif
         elseif(kw.eq.2.or.kw.eq.4)then
*           H-rich HG and CHeB temperature limit
            convradcomp = 10.d0**3.73d0
         elseif(kw.le.6)then
*           No limit or all other giant stars
*           (compare this large value to mass).
            convradcomp = 99999999.d0
         elseif(kw.eq.9)then
*           Mass limit for evolved He star
            convradcomp = 3.d0
         endif
*
*Magnetic Braking
*
         djmb = 0.d0
         if(htpmb.eq.0)then
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
                             djspint = djspint + djmb
            endif
         else
            if((menv.gt.0.d0.and.
     &            ((kw.le.1.and.mt.gt.0.35d0.and.
     &            mt.le.convradcomp).or.
     &            (kw.eq.2.and.teff.le.convradcomp).or.
     &            (kw.eq.4.and.teff.le.convradcomp).or.
     &            ((kw.eq.3).or.(kw.eq.5).or.
     &            (kw.eq.6)))))then
* MB given in Ivanova & Taam (2002)
               if(ospin.le.wx) djmb = kw3 * r**4.0d0 *
     &                                   (ospin/wsun)**3.0d0
               if(ospin.gt.wx) djmb = kw3 * r**4.0d0 *
     &                                   (ospin/wsun)**1.3d0 *
     &                                   (wx/wsun)**1.7d0
               djspint = djspint + djmb
            endif
         endif
         if(djmb.gt.tiny)then
            dtj = 0.03d0*jspin/ABS(djmb)
            dt = MIN(dt,dtj)
         endif
         if(kw.eq.13.and.pulsar.gt.0)then
*
* NS magnetic braking. PK.
* Single NS evolution.
*
            if(bdecayfac.eq.0)then
               if(B_0.eq.0.d0)then
                  B = 0.d0
               elseif((tphys-epoch-tacc).lt.tiny)then
                  B = B_0*EXP(-CK*bacc) + Bbot
               else
                  B = B_0*EXP(-(tphys-epoch-tacc)/bconst)*
     &                    EXP(-CK*bacc) + Bbot
               endif
            else
               if(B_0.eq.0.d0)then
                  B = 0.d0
               elseif((tphys-epoch-tacc).lt.tiny.and.
     &             bacc.eq.0.d0)then
                  B = B_0 + Bbot
               elseif((tphys-epoch-tacc).lt.tiny)then
                  B = B_0/(1.d0 + (bacc/1.0d-6)) + Bbot
               elseif(bacc.eq.0.d0)then
                  B = B_0*EXP(-(tphys-epoch-tacc)/bconst)
     &                    + Bbot
               else
                  B = B_0*EXP(-(tphys-epoch-tacc)/bconst)/
     &                (1.d0 + (bacc/1.0d-6)) + Bbot
               endif
            endif
            omdot = Kconst*B*B*ospin**3
            djmb = 0.4d0*mass*r*r*omdot
            djspint = djspint + djmb
* Consider update of time-stepping due to dj, i.e. dt = dj/(dj/dt).
            if(djmb.gt.tiny)then
               dtj = 0.1d0*(jspin/ABS(djmb))
               dt = MIN(dt,dtj)
            endif
         endif
         dtm = dt/1.0d+06

*
* Update mass and intrinsic spin (checking that the star is not spun
* past the equilibrium) and reset epoch for a MS (and possibly a HG) star.
*
*
         if(eqspin.gt.0.d0.and.ABS(dspint).gt.tiny)then
            if(intpol.eq.0)then
               if(dspint.ge.0.d0)then
                  dspint = MIN(dspint,(eqspin-ospin)/dt)
               else
                  dspint = MAX(dspint,(eqspin-ospin)/dt)
               endif
               djt = (k2str*(mt-mc)*r*r +
     &                k3*mc*rc*rc)*dspint
               djspint = djspint - djt
            endif
         endif
*
         jspin = MAX(1.0d-10,jspin - djspint*dt)
*
* Ensure that the star does not spin up beyond break-up.
*
         ospbru = twopi*SQRT(mt*aursun**3/r**3)
         jspbru = (k2str*(mt-mc)*r*r +
     &            k3*mc*rc*rc)*ospbru
         if((jspin.gt.jspbru.or.
     &      (jspin.eq.1.0d-10.and.intpol.gt.0)).and.
     &      ABS(dtm).gt.tiny)then !PDK add check for jspin intpol issues here.
* If rapidly spinning star (generally NS) pushed over the edge (spins up so
* much so that jspin would have become negative according to bse, so it sets
* it to 1.d0-10) in intpol then give it the stars maximum spin.
            mew = 1.d0
            if(djtx.gt.0.d0.and.jspin.gt.1.d0-10)then !dont go in here if jspin hit wall.
               mew = MIN(mew,(jspin - jspbru)/djtx)
            endif
             jspin = jspbru
         endif

*
* Update mass
*
         if(ABS(dms).gt.tiny)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tmold = tm
               tbgold = tscls(1)
               CALL star(kw,m0,mass,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch = tphys - epoch
                  endif
               else
                  epoch = tphys - ajhold*tm/tmold
               endif
            endif
* Update NS magnetic field owing to accretion, as a function of mass accreted. PK.
            if(kw.eq.13.and.pulsar.gt.0)then
               if(dms.lt.0.d0)then !negative dms is mass gained.
* When propeller ev. include .not.prop here...
                  b_mdot = dms/dt
                  if(b_mdot_lim.gt.0.d0.and.b_mdot.gt.b_mdot_lim)then
                     bacc = bacc - dms
                     tacc = tacc + dtm
                  elseif(b_mdot_lim.le.0.d0)then
                     bacc = bacc - dms
                     tacc = tacc + dtm
                  endif
               endif
            endif   
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*

         aj = tphys - epoch
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
         kwold = kw
         bhspin = 0.0
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r,lum,kw,mc,rc,menv,renv,k2,bhspin,1)
*
* If mass loss has occurred and no type change then check that we
* have indeed limited the radius change to 10%.
*
         if(kw.eq.kwold.and.dms.gt.0.d0)then
            mt2 = mt + dms
            dml = dms/dtm
            it = 0
 20         dr = r - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               it = it + 1
               if(it.eq.20.and.kw.eq.4) goto 30
               if(it.gt.30)then
                  WRITE(99,*)' DANGER1! ',it,kw,mass,dr,rm0
                  CALL FLUSH(99)
                  WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
                  CALL exit(0)
                  STOP
               endif
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(rm0,r)*dtdr
               if(it.ge.20) dtm = 0.5d0*dtm
               if(dtm.lt.1.0d-07*aj) goto 30
               dms = dtm*dml
               mt = mt2 - dms
               if(kw.le.2.or.kw.eq.7)then
                  mass = mt
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                  if(kw.eq.2)then
                     if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                        mass = m0
                     else
                        epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                               (tbgold - tmold)
                        epoch = tphys2 - epoch
                     endif
                  else
                     epoch = tphys2 - ajhold*tm/tmold
                  endif
               endif
               tphys = tphys2 + dtm
               aj = tphys - epoch
               mc = mc1
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                     r,lum,kw,mc,rc,menv,renv,k2,1)
               goto 20
            endif
 30         continue
         endif
*
* Initialize or adjust the spin of the star.
*
         if(j.eq.1)then
            if(tphys.lt.tiny.and.ospin.lt.0.001d0)then
               ospin = 45.35d0*vrotf(mt,0)/r
            endif
            jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         else
            jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
            ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         endif
*
* Test for changes in evolution type.
*
         if(j.eq.1.or.kw.ne.kwold)then
            if((kwold.le.12.and.(kw.eq.13.or.kw.eq.14)))then
               if(kw.eq.13.and.ecsn.gt.0.d0)then
                  if(kwold.le.6)then
                     if(mass.le.zpars(5))then
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                     endif
                  elseif(kwold.ge.7.and.kwold.le.9)then
                     if(mt.gt.ecsn_mlow.and.mt.le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                        if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                           sigma = -sigmahold/sigmadiv
                        elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                           sigma = sigmadiv
                        endif
                     endif
                  elseif(kwold.ge.10.or.kwold.eq.12)then
* AIC formation, will never happen here but...
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                  endif
               elseif(kw.eq.13.and.aic.gt.0)then
                  if(kwold.ge.10.or.kwold.eq.12)then
* AIC formation, will never happen here but...
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                  endif
               endif
               evolve_type = 14.d0 + FLOAT(k)
               teff = 1000.d0*((1130.d0*lum/
     &                       (r**2.d0))**(1.d0/4.d0))
               if(B_0.eq.0.d0)then !PK.
                  b01_bcm = 0.d0
               elseif(B_0.gt.0.d0.and.B.eq.0.d0)then
                  b01_bcm = B_0
               else
                  b01_bcm = B
               endif
               CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vk,k,
     &                   0.d0,fallback,sigmahold,kick_info,
     &                   disrupt)
            endif
*
* Force new NS or BH to have a birth spin peirod and magnetic field.
*
            if(kw.eq.13.or.kw.eq.14)then
               if(tphys-epoch.lt.tiny)then
                  if(kw.eq.13)then
*
 170                 u1 = ran3(idum1)
                     if(u1.ge.1.d0) goto 170
                     u2 = ran3(idum1)
                     s = sqrt(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
                     s = 0.7d0*s - 0.6d0
                     if(s.ge.0.013d0.or.s.le.-1.5d0) goto 170
                     ospin = (twopi*yearsc)/(10.d0**s)
 174                 u1 = ran3(idum1)
                     if(u1.ge.1.d0) goto 174
                     u2 = ran3(idum1)
                     s = sqrt(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
                     s = 0.68d0*s + 12.6d0
                     if(s.lt.11.5d0.or.s.gt.13.8d0) goto 174
                     B_0 = 10.d0**s
                     bacc = 0.d0 ! If it has been a NS before reset
                     tacc = 0.d0
                  endif
               else
                     ospin = 2.0d+08
               endif
               bacc = 0.d0 !make sure if its been a NS before its now a new one...
               tacc = 0.d0
               jspin = k3*rc*rc*mc*ospin
               sigma = sigmahold !reset sigma after possible ECSN kick dist.
            endif

*
* Record values for plotting and reset epoch.
*
            epoch = tphys - aj

            CALL WRITESPP(jp,tphys,evolve_type,
     &                    mass,kw,aj,tm,mc,r,m0,lum,
     &                    teff,rc,menv,renv,ospin,b01_bcm,
     &                    bacc,tacc,epoch,bhspin)
            CALL WRITESCM(ip,tphys,kw,m0,mass,lum,r,teff,
     &                    mc,rc,menv,renv,epoch,deltam,
     &                    ospin,b01_bcm,SN)
            jp = jp + 1
            ip = ip + 1
            if(kw.eq.15)then
               goto 90
            endif
         endif
*
         if(tphys.ge.tphysf)then
            CALL WRITESPP(jp,tphys,evolve_type,
     &                    mass,kw,aj,tm,mc,r,m0,lum,
     &                    teff,rc,menv,renv,ospin,B_0,
     &                    bacc,tacc,epoch,bhspin)
            goto 90
         endif
*
* Record radius and current age.
*
         rm0 = r
         ajhold = aj
         if(kw.ne.kwold) kwold = kw
         CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
*
* Check for type change.
*
         it = 0
         m0 = mass
         if((dtr-dtm).le.tol.and.kw.le.9)then
*
* Check final radius for too large a jump.
*
            aj = MAX(aj,aj*(1.d0-eps)+dtr)
            mc1 = mc
            CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r1,lum1,kw,mc1,rc1,menv1,renv1,k21,1)
            dr = r1 - rm0
            if(ABS(dr).gt.0.1d0*rm0)then
               dtm = dtr - ajhold*eps
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(r1,rm0)*dtdr
               goto 40
            else
               dtm = dtr
               goto 50
            endif
         endif
*
* Limit to a 10% increase in radius assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
 40      aj = ajhold + dtm
         mc1 = mc
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r1,lum1,kw,mc1,rc1,menv1,renv1,k21,1)
         dr = r1 - rm0
         it = it + 1
         if(it.eq.20.and.kw.eq.4) goto 50
         if(it.gt.30)then
            WRITE(99,*)' DANGER2! ',it,kw,mass,dr,rm0
            CALL FLUSH(99)
            WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
            CALL exit(0)
            STOP
         endif
         if(ABS(dr).gt.0.1d0*rm0)then
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm0,r1)*dtdr
            if(it.ge.20) dtm = 0.5d0*dtm
            goto 40
         endif
*
 50      continue
*
* Ensure that change of type has not occurred during radius check.
* This is rare but may occur for HG stars of ZAMS mass > 50 Msun.
*
         if(kw.ne.kwold)then
            kw = kwold
            mass = m0
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         endif

*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tsave-tphys)
*
 10   continue
*
 90   continue
*
      tphysf = tphys
      scm(ip+1,1) = -1.0
      spp(jp+1,1) = -1.0
      if(ip.ge.nv)then
         WRITE(99,*)' EVOLV1 ARRAY ERROR ',mass
         WRITE(*,*)' STOP: EVOLV1 ARRAY ERROR '
*         CALL exit(0)
*         STOP
      endif
*
      RETURN
      END
***
