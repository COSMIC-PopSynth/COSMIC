***
      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,ST_tide,
     &                  bhspin,kidx)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
*
*       H-R diagram for population I stars.
*       -----------------------------------
*
*       Computes the new mass, luminosity, radius & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout;
*       24th October 1995 to include metallicity;
*       14th November 1996 to include naked helium stars;
*       28th February 1997 to allow accretion induced supernovae.
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB
*
*
      integer kw,kwp,kidx
      INTEGER ST_tide
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 bhspin
      real*8 r,lum,mc,rc,menv,renv,k2
      real*8 mch,mlp,tiny
*      parameter(mch=1.44d0,mlp=12.d0,tiny=1.0d-14)
      parameter(mlp=12.d0,tiny=1.0d-14)
      real*8 mass0,mt0,mtc
      common /fall/fallback
      REAL*8 fallback
      REAL ran3
      EXTERNAL ran3
*
      real*8 mchold
*
      real*8 avar,bvar
      real*8 thook,thg,tbagb,tau,tloop,taul,tauh,tau1,tau2,dtau,texp
      real*8 lx,ly,dell,alpha,betahrdiag,eta
      real*8 rx,ry,delr,rzams,rtms,gammahrdiag,rmin,taumin,rg
      parameter(taumin=5.0d-08)
      real*8 mcmax,mcx,mcy,mcbagb,lambdahrdiag
      real*8 frac,kappa,sappa,alphap,polyfit
      real*8 am,xx,fac,rdgen,mew,lum0,kap,zeta,ahe,aco
      parameter(lum0=7.0d+04,kap=-0.5d0,ahe=4.d0,aco=16.d0)
*
      real*8 thookf,tblf
      real*8 lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      real*8 rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      real*8 rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      real*8 mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
      external thookf,tblf
      external lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      external rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      external rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      external mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old; output: new value).
*       AJ      Current age in Myr.
*       MT      Current mass in solar units (used for R).
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       TSCLS   Time scale for different stages.
*       LUMS    Characteristic luminosity.
*       GB      Giant Branch parameters
*       ZPARS   Parameters for distinguishing various mass intervals.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (0 - 15).
*       MC      Core mass.
*       ---------------------------------------------------------------------
*
*
* Make evolutionary changes to stars that have not reached KW > 5.
*
      mch = 1.44d0 !set here owing to AIC ECSN model.
*
      mass0 = mass
C      if(mass0.gt.100.d0) mass = 100.d0
      mt0 = mt
C      if(mt0.gt.100.d0) mt = 100.d0
*
      if(kw.gt.6) goto 90
*
      tbagb = tscls(2) + tscls(3)
      thg = tscls(1) - tm
*
      rzams = rzamsf(mass)
      rtms = rtmsf(mass)
*
      if(aj.lt.tscls(1))then
*
*        Either on MS or HG
*
         rg = rgbf(mt,lums(3))
*
         if(aj.lt.tm)then
*
*           Main sequence star.
*
            mc = 0.d0
            tau = aj/tm
            thook = thookf(mass)*tscls(1)
            zeta = 0.01d0
            tau1 = MIN(1.d0,aj/thook)
            tau2 = MAX(0.d0,
     &             MIN(1.d0,(aj-(1.d0-zeta)*thook)/(zeta*thook)))
*
            dell = lhookf(mass,zpars(1))
            dtau = tau1**2 - tau2**2
            alpha = lalphf(mass)
            betahrdiag = lbetaf(mass)
            eta = lnetaf(mass)
            lx = LOG10(lums(2)/lums(1))
            if(tau.gt.taumin)then
               xx = alpha*tau + betahrdiag*tau**eta +
     &              (lx - alpha - betahrdiag)*tau**2 - dell*dtau
            else
               xx = alpha*tau + (lx - alpha)*tau**2 - dell*dtau
            endif
            lum = lums(1)*10.d0**xx
*
            delr = rhookf(mass,zpars(1))
            dtau = tau1**3 - tau2**3
            alpha = ralphf(mass)
            betahrdiag = rbetaf(mass)
            gammahrdiag = rgammf(mass)
            rx = LOG10(rtms/rzams)
* Note that the use of taumin is a slightly pedantic attempt to
* avoid floating point underflow. It IS overkill!
            if(tau.gt.taumin)then
               xx = alpha*tau + betahrdiag*tau**10 +
     &              gammahrdiag*tau**40 + (rx - alpha - betahrdiag -
     &              gammahrdiag)*tau**3 - delr*dtau
            else
               xx = alpha*tau + (rx - alpha)*tau**3 - delr*dtau
            endif
            r = rzams*10.d0**xx
*
            if(mass.lt.(zpars(1)-0.3d0))then
               kw = 0
* This following is given by Chris for low mass MS stars which will be
* substantially degenerate. We need the Hydrogen abundance, X, which we
* calculate from Z assuming that the helium abundance, Y, is calculated
* according to Y = 0.24 + 2*Z
               rdgen = 0.0258d0*((1.d0+zpars(11))**(5.d0/3.d0))*
     &                          (mass**(-1.d0/3.d0))
               r = MAX(rdgen,r)
            else
               kw = 1
            endif
* planets
            if(mass.lt.0.005d0.and.mass.ge.tiny)then
               r = 0.16d0
            endif
*
         else
*
*           Star is on the HG
*
            mcx = mc
            if(mass.le.zpars(2))then
               mc = mcgbf(lums(3),GB,lums(6))
            elseif(mass.le.zpars(3))then
               mc = mcheif(mass,zpars(2),zpars(9))
            else
               mc = mcheif(mass,zpars(2),zpars(10))
            endif
            eta = mctmsf(mass)
            tau = (aj - tm)/thg
            mc = ((1.d0 - tau)*eta + tau)*mc
            mc = MAX(mc,mcx)
*
* Test whether core mass has reached total mass.
*
            if(mc.ge.mt)then
               aj = 0.d0
               if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
                  mc = 0.d0
                  mass = mt
                  kw = 7
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               else
*
* Zero-age helium white dwarf.
*
                  mc = mt
                  mass = mt
                  kw = 10
               endif
            else
               lum = lums(2)*(lums(3)/lums(2))**tau
               if(mass.le.zpars(3))then
                  rx = rg
               else
* He-ignition and end of HG occur at Rmin
                  rmin = rminf(mass)
                  ry = ragbf(mt,lums(4),zpars(2))
                  rx = MIN(rmin,ry)
                  if(mass.le.mlp)then
                     texp = log(mass/mlp)/log(zpars(3)/mlp)
                     rx = rg
                     rx = rmin*(rx/rmin)**texp
                  endif
                  tau2 = tblf(mass,zpars(2),zpars(3))
                  if(tau2.lt.tiny) rx = ry
               endif
               r = rtms*(rx/rtms)**tau
               kw = 2
            endif
*
         endif
*
* Now the GB, CHeB and AGB evolution.
*
      elseif(aj.lt.tscls(2))then
*
*        Red Giant.
*
         kw = 3
         lum = lgbtf(aj,GB(1),GB,tscls(4),tscls(5),tscls(6))
         if(mass.le.zpars(2))then
* Star has a degenerate He core which grows on the GB
            mc = mcgbf(lum,GB,lums(6))
         else
* Star has a non-degenerate He core which may grow, but
* only slightly, on the GB
            tau = (aj - tscls(1))/(tscls(2) - tscls(1))
            mcx = mcheif(mass,zpars(2),zpars(9))
            mcy = mcheif(mass,zpars(2),zpars(10))
            mc = mcx + (mcy - mcx)*tau
         endif
         r = rgbf(mt,lum)
         rg = r
         if(mc.ge.mt)then
            aj = 0.d0
            if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
               mc = 0.d0
               mass = mt
               kw = 7
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            else
*
* Zero-age helium white dwarf.
*
               mc = mt
               mass = mt
               kw = 10
            endif
         endif
*
      elseif(aj.lt.tbagb)then
*
*       Core helium burning star.
*
         if(kw.eq.3.and.mass.le.zpars(2))then
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2)
         endif
         if(mass.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(mass,zpars(2),zpars(10))
         endif
         tau = (aj - tscls(2))/tscls(3)
         mc = mcx + (mcagbf(mass) - mcx)*tau
*
         if(mass.le.zpars(2))then
            lx = lums(5)
            ly = lums(7)
            rx = rzahbf(mt,mc,zpars(2))
            rg = rgbf(mt,lx)
            rmin = rg*zpars(13)**(mass/zpars(2))
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            ry = ragbf(mt,ly,zpars(2))
            if(rmin.lt.rx)then
               taul = (log(rx/rmin))**(1.d0/3.d0)
            else
               rmin = rx
               taul = 0.d0
            endif
            tauh = (log(ry/rmin))**(1.d0/3.d0)
            tau2 = taul*(tau - 1.d0) + tauh*tau
            r = rmin*exp(abs(tau2)**3)
            rg = rg + tau*(ry - rg)
            lum = lx*(ly/lx)**(tau**texp)
         elseif(mass.gt.zpars(3))then
*
* For HM stars He-ignition takes place at Rmin in the HG, and CHeB
* consists of a blue phase (before tloop) and a RG phase (after tloop).
*
            tau2 = tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            rmin = rminf(mass)
            rg = rgbf(mt,lums(4))
            rx = ragbf(mt,lums(4),zpars(2))
            rmin = MIN(rmin, rx)
            if(mass.le.mlp) then
               texp = log(mass/mlp)/log(zpars(3)/mlp)
               rx = rg
               rx = rmin*(rx/rmin)**texp
            else
               rx = rmin
            endif
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            lum = lums(4)*(lums(7)/lums(4))**(tau**texp)
            if(aj.lt.tloop)then
               ly = lums(4)*(lums(7)/lums(4))**(tau2**texp)
               ry = ragbf(mt,ly,zpars(2))
               taul = 0.d0
               if(ABS(rmin-rx).gt.tiny)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               endif
               tauh = 0.d0
               if(ry.gt.rmin) tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tscls(2))/(tau2*tscls(3))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rg + tau*(ry - rg)
            else
               r = ragbf(mt,lum,zpars(2))
               rg = r
            endif
         else
*
* For IM stars CHeB consists of a RG phase (before tloop) and a blue
* loop (after tloop).
*
            tau2 = 1.d0 - tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            if(aj.lt.tloop)then
               tau = (tloop - aj)/(tau2*tscls(3))
               lum = lums(5)*(lums(4)/lums(5))**(tau**3)
               r = rgbf(mt,lum)
               rg = r
            else
               lx = lums(5)
               ly = lums(7)
               rx = rgbf(mt,lx)
               rmin = rminf(mt)
               texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
               ry = ragbf(mt,ly,zpars(2))
               if(rmin.lt.rx)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               else
                  rmin = rx
                  taul = 0.d0
               endif
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tloop)/(tscls(3) - (tloop - tscls(2)))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rx + tau*(ry - rx)
               lum = lx*(ly/lx)**(tau**texp)
            endif
         endif
*
* Test whether core mass exceeds total mass.
*
         if(mc.ge.mt)then
*
* Evolved MS naked helium star.
*
            kw = 7
            xx = (aj - tscls(2))/tscls(3)
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = xx*tm
         else
            kw = 4
         endif
*
      else
*
*        Asymptotic Red Giant.
*
* On the AGB the He core mass remains constant until at Ltp it
* is caught by the C core mass and they grow together.
*
         mcbagb = mcagbf(mass)
         mcx = mcgbtf(tbagb,GB(8),GB,tscls(7),tscls(8),tscls(9))
         mcmax = MAX(MAX(mch,0.773d0*mcbagb-0.35d0),1.05d0*mcx)
*
         if(aj.lt.tscls(13))then
            mcx = mcgbtf(aj,GB(8),GB,tscls(7),tscls(8),tscls(9))
            mc = mcbagb
            lum = lmcgbf(mcx,GB)
            if(mt.le.mc)then
*
* Evolved naked helium star as the envelope is lost but the
* star has not completed its interior burning. The star becomes
* a post-HeMS star.
*
               kw = 9
               mt = mc
               mass = mt
               mc = mcx
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(mc.le.GB(7))then
                  aj = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                            (mc**(1.d0-GB(5)))
               else
                  aj = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                            (mc**(1.d0-GB(6)))
               endif
               aj = MAX(aj,tm)
               goto 90
            else
               kw = 5
            endif
         else
            kw = 6
            mc = mcgbtf(aj,GB(2),GB,tscls(10),tscls(11),tscls(12))
            lum = lmcgbf(mc,GB)
*
* Approximate 3rd Dredge-up on AGB by limiting Mc.
*
            lambdahrdiag = MIN(0.9d0,0.3d0+0.001d0*mass**5)
            tau = tscls(13)
            mcx = mcgbtf(tau,GB(2),GB,tscls(10),tscls(11),tscls(12))
            mcy = mc
            mc = mc - lambdahrdiag*(mcy-mcx)
            mcx = mc
            mcmax = MIN(mt,mcmax)
         endif
         r = ragbf(mt,lum,zpars(2))
         rg = r
*
* Mc,x represents the C core mass and we now test whether it
* exceeds either the total mass or the maximum allowed core mass.
*
         if(mcmax-mcx.lt.tiny)then
            aj = 0.d0
            mc = mcmax
            if(mc.lt.mch)then
               if(ifflag.ge.1)then
*
* Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800.
*
                  if(zpars(14).ge.1.0d-08)then
                     mc = MIN(0.36d0+0.104d0*mass,0.58d0+0.061d0*mass)
                     mc = MAX(0.54d0+0.042d0*mass,mc)
                     if(mass.lt.1.d0) mc = 0.46d0
                  else
                     mc = MIN(0.29d0+0.178d0*mass,0.65d0+0.062d0*mass)
                     mc = MAX(0.54d0+0.073d0*mass,mc)
                  endif
                  mc = MIN(mch,mc)
               endif
*
               mt = mc
               if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
                  kw = 11
               elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then !double check what this should be. should be ecsn_mlow. Remember need to add option if ecsn = 0 (i.e. no ECSN!!!)
*
* Zero-age Carbon/Oxygen White Dwarf
*
                  kw = 11
               elseif(ecsn.gt.0.d0.and.mcbagb.ge.ecsn_mlow.and.
     &                mcbagb.le.ecsn.and.mc.lt.1.08d0)then
                  kw = 11
*               elseif(mcbagb.ge.1.6d0.and.mcbagb.le.2.5d0.and.
*                      mc.lt.1.08d0)then !can introduce this into code at some point.
*                  kw = 11

               else
*
* Zero-age Oxygen/Neon White Dwarf
*
                  kw = 12
               endif
               mass = mt
*
            else
               if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  kw = 15
                  aj = 0.d0
                  mt = 0.d0
                  lum = 1.0d-10
                  r = 1.0d-10
               elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  kw = 15
                  aj = 0.d0
                  mt = 0.d0
                  lum = 1.0d-10
                  r = 1.0d-10
               else
                  if(nsflag.eq.0)then
                     mt = 1.17d0 + 0.09d0*mc
                  elseif(nsflag.eq.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407.
*
                     if(mc.lt.2.5d0)then
                        mcx = 0.161767d0*mc + 1.067055d0
                     else
                        mcx = 0.314154d0*mc + 0.686088d0
                     endif
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                        fallback = (mc - 5.d0)/2.6d0
                     elseif(mc.gt.7.60)then
                        fallback = 1.d0
                     endif
                  elseif(nsflag.eq.2)then
*
* Use NS/BH masses given by Belczynski+08. PK.
*
                     !First calculate the proto-core mass
                     if(ecsn.gt.0.d0.and.mcbagb.le.ecsn)then
                        mcx = 1.38d0
                     elseif(ecsn.eq.0.d0.and.mcbagb.le.2.25d0)then !this should be ecsn, unless ecsn=0
*                     if(mcbagb.le.2.35d0)then
                        mcx = 1.38d0
*                     elseif(mc.lt.4.29d0)then
                     elseif(mc.lt.4.82d0)then
                        mcx = 1.5d0
*                     elseif(mc.ge.4.29d0.and.mc.lt.6.31d0)then
                     elseif(mc.ge.4.82d0.and.mc.lt.6.31d0)then
                        mcx = 2.11d0
                     elseif(mc.ge.6.31d0.and.mc.lt.6.75d0)then
                        mcx = 0.69*mc - 2.26d0
                     elseif(mc.ge.6.75d0)then
                        mcx = 0.37*mc - 0.07d0
                     endif
                     !now calculate the remnant mass after fallback
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (0.378d0*mc - 1.889d0)*(mt - mcx)
                        fallback = (0.378d0*mc - 1.889d0)
                     elseif(mc.gt.7.60)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(nsflag.eq.3)then
*
* Use the "Rapid" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    For this, we just set the proto-core mass to one
                     mcx = 1.d0
                     if(ecsn.gt.0.d0.and.mcbagb.le.ecsn)then
                        mcx = 1.38d0
                     elseif(ecsn.eq.0.d0.and.mcbagb.le.2.25d0)then !this should be ecsn, unless ecsn=0
                        mcx = 1.38d0
                     endif
                     if(mc.le.2.5d0)then
                        mt = mcx + 0.2d0
                        fallback = 0.d0
                     elseif(mc.le.6.d0)then
                        fallback = (0.286d0*mc - 0.514d0) / (mt - mcx)
                        mt = mcx + 0.286d0*mc - 0.514d0
                     elseif(mc.le.7.d0)then
                        fallback = 1.d0
                     elseif(mc.le.11.d0)then
                        avar = 0.25d0 - (1.275 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.gt.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(nsflag.eq.4)then
*
* Use the "Delayed" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    For this, we just set the proto-core mass to one
                     if(mc.le.3.5d0)then
                        mcx = 1.2d0
                     elseif(mc.le.6.d0)then
                        mcx = 1.3d0
                     elseif(mc.le.11.d0)then
                        mcx = 1.4d0
                     elseif(mc.gt.11.d0)then
                        mcx = 1.6d0
                     endif
                     if(mc.lt.2.5d0)then
                        mt = mcx + 0.2
                        fallback = 0.d0
                     elseif(mc.lt.3.5d0)then
                        fallback = (0.5d0 * mc - 1.05d0) / (mt - mcx)
                        mt = mcx + 0.5d0 * mc - 1.05d0
                     elseif(mc.lt.11.d0)then
                        avar = 0.133d0 - (0.093d0 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.ge.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  endif

                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
* Convert baryonic mass to gravitational mass (Lattimer & Yahil 1989)
                     if(nsflag.ge.2) mt = 2.108167d0*SQRT(10.d0+3.d0*mt)
     &                                    - 6.6666667d0
                  else
*
* Zero-age Black hole
*
                     kw = 14

* CLR - (Pulsational) Pair-Instability Supernova

* Belczynski+2016 prescription: just shrink any BH with a He core mass
* between 45 and 65 solar masses (provided the pisn flag is set at 45),
* and blow up anything between 65 and 135 solar masses.  
* Cheap, but effective
                     if(pisn.gt.0)then
                        if(mcbagb.ge.pisn.and.mcbagb.lt.65.d0)then
                           mt = pisn
                           mc = pisn
                           pisn_track(kidx)=6
                        elseif(mcbagb.ge.65.d0.and.mcbagb.lt.135.d0)then
                           mt = 0.d0
                           mc = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif
* The Spera+Mapelli2017 prescription is a tad more sophisticated:
* complex fitting formula to Stan Woosley's PSN models.  HOWEVER, these
* were done using the ZAMS mass/core mass/remnant mass relationships for
* SEVN, not BSE.  In other words, I woud be careful using this (and in
* practice, it doesn't vary that much from Belczynski's prescription,
* since the He core masses are the same in both)
                     elseif(pisn.eq.-1)then
                        frac = mcbagb/mt
                        kappa = 0.67d0*frac + 0.1d0
                        sappa = 0.5226d0*frac - 0.52974d0
                        if(mcbagb.le.32.d0)then
                           alphap = 1.0d0
                        elseif(frac.lt.0.9d0.and.mcbagb.le.37.d0)then
                           alphap = 0.2d0*(kappa-1.d0)*mcbagb +
     &                              0.2d0*(37.d0 - 32.d0*kappa)
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mcbagb.le.60.d0)then
                           alphap = kappa
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mcbagb.lt.64.d0)then
                           alphap = kappa*(-0.25d0)*mcbagb + kappa*16.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.37.d0)then
                           alphap = sappa*(mcbagb - 32.d0) + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.56.d0.and.
     &                         sappa.lt.-0.034168d0)then
                           alphap = 5.d0*sappa + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.56.d0.and.
     &                         sappa.ge.-0.034168d0)then
                           alphap = (-0.1381d0*frac + 0.1309d0)*
     &                              (mcbagb - 56.d0) + 0.82916d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.lt.64.d0)then
                           alphap = -0.103645d0*mcbagb + 6.63328d0
                           pisn_track(kidx)=6
                        elseif(mcbagb.ge.64.d0.and.mcbagb.lt.135.d0)then
                           alphap = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        elseif(mcbagb.ge.135.d0)then
                           alphap = 1.0d0
                        endif
                        mt = alphap*mt

* Fit (8th order polynomial) to Table 1 in Marchant+2018.
                     elseif(pisn.eq.-2)then
                        if(mcbagb.ge.27.69d0.and.mcbagb.le.54.48d0)then
                           polyfit = -4.30343374d5
     &                            + 9.02795937d4*mcbagb
     &                            - 8.22480314d3*mcbagb**2d0
     &                            + 4.25048530d2*mcbagb**3d0
     &                            - 1.36291200d1*mcbagb**4d0
     &                            + 2.77684136d-1*mcbagb**5d0
     &                            - 3.51104067d-3*mcbagb**6d0
     &                            + 2.51918414d-5*mcbagb**7d0
     &                            - 7.85424404d-8*mcbagb**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mcbagb.gt.54.48d0.and.
     &                         mcbagb.lt.113.29d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif

* Fit (8th order polynomial) to Table 5 in Woosley2019.
                     elseif(pisn.eq.-3)then
                        if(mcbagb.ge.29.53d0.and.mcbagb.le.60.12d0)then
                           polyfit = -3.14610870d5
     &                            + 6.13699616d4*mcbagb
     &                            - 5.19249710d3*mcbagb**2d0
     &                            + 2.48914888d2*mcbagb**3d0
     &                            - 7.39487537d0*mcbagb**4d0
     &                            + 1.39439936d-1*mcbagb**5d0
     &                            - 1.63012111d-3*mcbagb**6d0
     &                            + 1.08052344d-5*mcbagb**7d0
     &                            - 3.11019088d-8*mcbagb**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mcbagb.gt.60.12d0.and.
     &                         mcbagb.lt.135.d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif
                     endif

* Convert baryonic mass to gravitational mass (approx for BHs)
                     if(nsflag.ge.2) mt = 0.9d0*mt
                  endif
               endif
            endif
         endif
*
      endif
*
 90   continue
*
      if(kw.ge.7.and.kw.le.9)then
*
* Naked Helium Star
*
         rzams = rzhef(mt)
         rx = rzams
         if(aj.lt.tm)then
*
* Main Sequence
*
            kw = 7
            tau = aj/tm
            am = MAX(0.d0,0.85d0-0.08d0*mass)
            lum = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
            am = MAX(0.d0,0.4d0-0.22d0*LOG10(mt))
            r = rx*(1.d0+am*(tau-tau**6))
            rg = rx
* Star has no core mass and hence no memory of its past
* which is why we subject mass and mt to mass loss for
* this phase.
            mc = 0.d0
            if(mt.lt.zpars(10)) kw = 10
         else
*
* Helium Shell Burning
*
            kw = 8
            lum = lgbtf(aj,GB(8),GB,tscls(4),tscls(5),tscls(6))
            r = rhehgf(mt,lum,rx,lums(2))
            rg = rhegbf(lum)
            if(r.ge.rg)then
               kw = 9
               r = rg
            endif
            mc = mcgbf(lum,GB,lums(6))
            mtc = MIN(mt,1.45d0*mt-0.31d0)
            mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
            if(mcmax-mc.lt.tiny)then
               aj = 0.d0
               mc = mcmax
               if(mc.lt.mch)then
                  if(ecsn.gt.0.d0.and.mass.lt.ecsn_mlow)then
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  elseif(ecsn.eq.0.d0.and.mass.lt.1.6d0)then
*
* Zero-age Carbon/Oxygen White Dwarf
*
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  elseif(ecsn.gt.0.d0.and.mass.gt.ecsn_mlow.and.
     &                   mass.le.ecsn.and.mc.le.1.08d0)then
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  else
*
* Zero-age Oxygen/Neon White Dwarf
*
                     mt = mc
                     kw = 12
                  endif
                  mass = mt
               else
                  if(ecsn.gt.0.d0.and.mass.lt.ecsn_mlow)then
                     kw = 15
                     aj = 0.d0
                     mt = 0.d0
                     lum = 1.0d-10
                     r = 1.0d-10
                  elseif(ecsn.eq.0.d0.and.mass.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                     kw = 15
                     aj = 0.d0
                     mt = 0.d0
                     lum = 1.0d-10
                     r = 1.0d-10
                  else
                     if(nsflag.eq.0)then
                        mt = 1.17d0 + 0.09d0*mc
                     elseif(nsflag.eq.1)then
                        if(mc.lt.2.5d0)then
                           mcx = 0.161767d0*mc + 1.067055d0
                        else
                           mcx = 0.314154d0*mc + 0.686088d0
                        endif
                        if(mc.le.5.d0)then
                           mt = mcx
                           fallback = 0.d0
                        elseif(mc.lt.7.6d0)then
                           mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                           fallback = (mc - 5.d0)/2.6d0
                        endif
                        if(mc.gt.7.60) fallback = 1.d0
                     elseif(nsflag.eq.2)then
*
* Use NS/BH masses given by Belczynski+08. PK.
*
                     !First calculate the proto-core mass
                     if(ecsn.gt.0.d0.and.mc.le.ecsn)then
                        mcx = 1.38d0
                     elseif(ecsn.eq.0.d0.and.mc.le.2.25d0)then !this should be ecsn, unless ecsn=0
                        mcx = 1.38d0
                     elseif(mc.lt.4.82d0)then
                        mcx = 1.5d0
                     elseif(mc.ge.4.82d0.and.mc.lt.6.31d0)then
                        mcx = 2.11d0
                     elseif(mc.ge.6.31d0.and.mc.lt.6.75d0)then
                        mcx = 0.69*mc - 2.26d0
                     elseif(mc.ge.6.75d0)then
                        mcx = 0.37*mc - 0.07d0
                     endif
                     !now calculate the remnant mass after fallback
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (0.378d0*mc - 1.889d0)*(mt - mcx)
                        fallback = (0.378d0*mc - 1.889d0)
                     elseif(mc.gt.7.60)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(nsflag.eq.3)then
*
* Use the "Rapid" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    For this, we just set the proto-core mass to one
                     mcx = 1.d0
                     if(ecsn.gt.0.d0.and.mc.le.ecsn)then
                        mcx = 1.38d0
                     elseif(ecsn.eq.0.d0.and.mc.le.2.25d0)then !this should be ecsn, unless ecsn=0
                        mcx = 1.38d0
                     endif
                     if(mc.le.2.5d0)then
                        mt = mcx + 0.2d0
                        fallback = 0.d0
                     elseif(mc.le.6.d0)then
                        fallback = (0.286d0*mc - 0.514d0) / (mt - mcx)
                        mt = mcx + 0.286d0*mc - 0.514d0
                     elseif(mc.le.7.d0)then
                        fallback = 1.d0
                     elseif(mc.le.11.d0)then
                        avar = 0.25d0 - (1.275 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.gt.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(nsflag.eq.4)then
*
* Use the "Delayed" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    For this, we just set the proto-core mass to one
                     if(mc.le.3.5d0)then
                        mcx = 1.2d0
                     elseif(mc.le.6.d0)then
                        mcx = 1.3d0
                     elseif(mc.le.11.d0)then
                        mcx = 1.4d0
                     elseif(mc.gt.11.d0)then
                        mcx = 1.6d0
                     endif
                     if(mc.lt.2.5d0)then
                        mt = mcx + 0.2
                        fallback = 0.d0
                     elseif(mc.lt.3.5d0)then
                        fallback = (0.5d0 * mc - 1.05d0) / (mt - mcx)
                        mt = mcx + 0.5d0 * mc - 1.05d0
                     elseif(mc.lt.11.d0)then
                        avar = 0.133d0 - (0.093d0 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.ge.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  endif

                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
* Convert baryonic mass to gravitational mass (Lattimer & Yahil 1989)
                     if(nsflag.ge.2) mt = 2.108167d0*SQRT(10.d0+3.d0*mt)
     &                                    - 6.6666667d0
                  else
*
* Zero-age Black hole
*
                     kw = 14

* CLR - (Pulsational) Pair-Instability Supernova

* Belczynski+2016 prescription: just shrink any BH with a He core mass
* between 45 and 65 solar masses, and blow up anything between 65 and
* 135 solar masses.  Cheap, but effective
                     if(pisn.gt.0)then
                        if(mc.ge.pisn.and.mc.lt.65.d0)then
                           mt = pisn
                           mc = pisn
                           pisn_track(kidx)=6
                        elseif(mc.ge.65.d0.and.mc.lt.135.d0)then
                           mt = 0.d0
                           mc = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif
* The Spera+Mapelli2017 prescription is a tad more sophisticated:
* complex fitting formula to Stan Woosley's PSN models.  HOWEVER, these
* were done using the ZAMS mass/core mass/remnant mass relationships for
* SEVN, not BSE.  In other words, I woud be careful using this (and in
* practice, it doesn't vary that much from Belczynski's prescription,
* since the He core masses are the same in both)
* Mario said this prescription works here as well.
                     elseif(pisn.eq.-1)then
                        frac = mc/mt
                        kappa = 0.67d0*frac + 0.1d0
                        sappa = 0.5226d0*frac - 0.52974d0
                        if(mc.le.32.d0)then
                           alphap = 1.0d0
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mc.le.37.d0)then
                           alphap = 0.2d0*(kappa-1.d0)*mc +
     &                              0.2d0*(37.d0 - 32.d0*kappa)
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mc.le.60.d0)then
                           alphap = kappa
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mc.lt.64.d0)then
                           alphap = kappa*(-0.25d0)*mc+ kappa*16.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mc.le.37.d0)then
                           alphap = sappa*(mc- 32.d0) + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mc.le.56.d0.and.
     &                         sappa.lt.-0.034168d0)then
                           alphap = 5.d0*sappa + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mc.le.56.d0.and.
     &                         sappa.ge.-0.034168d0)then
                           alphap = (-0.1381d0*frac + 0.1309d0)*
     &                              (mc- 56.d0) + 0.82916d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mc.lt.64.d0)then
                           alphap = -0.103645d0*mc+ 6.63328d0
                           pisn_track(kidx)=6
                        elseif(mc.ge.64.d0.and.mc.lt.135.d0)then
                           alphap = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        elseif(mc.ge.135.d0)then
                           alphap = 1.0d0
                        endif
                        mt = alphap*mt


* Fit (8th order polynomial) to Table 1 in Marchant+2018.
                     elseif(pisn.eq.-2)then
                        if(mc.ge.27.69d0.and.mc.le.54.48d0)then
                           polyfit = -4.30343374d5
     &                            + 9.02795937d4*mc
     &                            - 8.22480314d3*mc**2d0
     &                            + 4.25048530d2*mc**3d0
     &                            - 1.36291200d1*mc**4d0
     &                            + 2.77684136d-1*mc**5d0
     &                            - 3.51104067d-3*mc**6d0
     &                            + 2.51918414d-5*mc**7d0
     &                            - 7.85424404d-8*mc**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mc.gt.54.48d0.and.mc.lt.113.29d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif

* Fit (8th order polynomial) to Table 5 in Woosley2019.
                     elseif(pisn.eq.-3)then
                        if(mc.ge.29.53d0.and.mc.le.60.12d0)then
                           polyfit = -3.14610870d5
     &                            + 6.13699616d4*mc
     &                            - 5.19249710d3*mc**2d0
     &                            + 2.48914888d2*mc**3d0
     &                            - 7.39487537d0*mc**4d0
     &                            + 1.39439936d-1*mc**5d0
     &                            - 1.63012111d-3*mc**6d0
     &                            + 1.08052344d-5*mc**7d0
     &                            - 3.11019088d-8*mc**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mc.gt.60.12d0.and.mc.lt.135.d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif
                     endif


* Convert baryonic mass to gravitational mass (approx for BHs)
                     if(nsflag.ge.2) mt = 0.9d0*mt
                     endif
                  endif
               endif
            endif
         endif
      endif
*
      if(kw.ge.10.and.kw.le.12)then
*
*        White dwarf.
*
         mc = mt
         mchold = mch
         if(ecsn.gt.0.d0.and.kw.eq.12) mch = 1.38d0
         if(mc.ge.mch)then
*
* Accretion induced supernova with no remnant
* unless WD is ONe in which case we assume a NS
* of minimum mass is the remnant.
*
            if(kw.eq.12)then
               kw = 13
               aj = 0.d0
               mt = 1.3d0
               if(ecsn.gt.0.d0)then
                  mt = 1.38d0
                  mt = 0.9d0*mt !in ST this is a quadratic, will add in later.
               endif
            else
               kw = 15
               aj = 0.d0
               mt = 0.d0
               lum = 1.0d-10
               r = 1.0d-10
            endif
         else
*
            if(kw.eq.10)then
               xx = ahe
            else
               xx = aco
            endif
*
            if(wdflag.eq.0)then
*
* Mestel cooling
*
               lum = 635.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.4d0
*
            elseif(wdflag.ge.1)then
*
* modified-Mestel cooling
*
               if(aj.lt.9000.0)then
                  lum = 300.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.18d0
               else
                  fac = (9000.1d0*xx)**5.3d0
                  lum = 300.d0*fac*mt*zpars(14)/(xx*(aj+0.1d0))**6.48d0
               endif
*
            endif
*
            r = 0.0115d0*SQRT(MAX(1.48204d-06,(mch/mt)**(2.d0/3.d0)
     &                                      - (mt/mch)**(2.d0/3.d0)))
            r = MIN(0.1d0,r)
            if(mt.lt.0.0005d0) r = 0.09d0
            if(mt.lt.0.000005d0) r = 0.009d0
*
         endif
         mch = mchold !added for AIC ECSN stuff.
      endif
*
      if(kw.eq.13)then
*
*        Neutron Star.
*
         mc = mt
         if(mc.gt.mxns)then
*
* Accretion induced Black Hole?
*
            kw = 14
            aj = 0.d0
         else
            lum = 0.02d0*(mt**0.67d0)/(MAX(aj,0.1d0))**2
            r = 1.4d-05
         endif
      endif
*
      if(kw.eq.14)then
*
*        Black hole
*
         mc = mt
         lum = 1.0d-10
         r = 4.24d-06*mt
      endif
*
* Calculate the core radius and the luminosity and radius of the
* remnant that the star will become.
*
      tau = 0.d0
      if(kw.le.1.or.kw.eq.7)then
         rc = 0.d0
      elseif(kw.le.3)then
         if(mass.gt.zpars(2))then
            lx = lzhef(mc)
            rx = rzhef(mc)
            rc = rx
         else
            if(wdflag.eq.0)then
               lx = 635.d0*mc*zpars(14)/((ahe*0.1d0)**1.4d0)
            elseif(wdflag.ge.1)then
               lx = 300.d0*mc*zpars(14)/((ahe*0.1d0)**1.18d0)
            endif
            rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &           (mch/mc)**(2.d0/3.d0)-(mc/mch)**(2.d0/3.d0)))
            rc = 5.d0*rx
         endif
      elseif(kw.eq.4)then
         tau = (aj - tscls(2))/tscls(3)
         kwp = 7
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         am = MAX(0.d0,0.85d0-0.08d0*mc)
         lx = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
         rx = rzhef(mc)
         am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
         rx = rx*(1.d0+am*(tau-tau**6))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.eq.5)then
         kwp = 9
         if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         lx = lmcgbf(mcx,GB)
         if(tau.lt.1.d0) lx = lums(2)*(lx/lums(2))**tau
         rx = rzhef(mc)
         rx = MIN(rhehgf(mc,lx,rx,lums(2)),rhegbf(lx))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.le.9)then
         if(wdflag.eq.0)then
            lx = 635.d0*mc*zpars(14)/((aco*0.1d0)**1.4d0)
         elseif(wdflag.ge.1)then
            lx = 300.d0*mc*zpars(14)/((aco*0.1d0)**1.18d0)
         endif
         rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &        (mch/mc)**(2.d0/3.d0) - (mc/mch)**(2.d0/3.d0)))
         rc = 5.d0*rx
      else
         rc = r
         menv = 1.0d-10
         renv = 1.0d-10
         k2 = 0.21d0
      endif
*
* Perturb the luminosity and radius due to small envelope mass.
*
      if(kw.ge.2.and.kw.le.9.and.kw.ne.7)then
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
         if(kw.ge.8) mew = ((mtc-mc)/mtc)*5.d0
         if(mew.lt.1.d0)then
            xx = lpertf(mt,mew)
            lum = lx*(lum/lx)**xx
            if(r.le.rx)then
               xx = 0.d0
            else
               xx = rpertf(mt,mew,r,rx)
            endif
            r = rx*(r/rx)**xx
         endif
         rc = MIN(rc,r)
      endif
*
* Calculate mass and radius of convective envelope, and envelope
* gyration radius.
*
      if(kw.lt.10)then
         CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),
     &              lums(4),rzams,rtms,rg,menv,renv,k2)
      endif
*
      if(ST_tide.gt.0)then
         if(kw.le.2.or.kw.eq.7.or.kw.ge.10)then
            if(mt.le.1.d0)then
               k2 = 0.205d0
            else
               k2 = 0.075d0
            endif
         else
           k2 = 0.1d0
         endif
      endif
*
C      if(mass.gt.99.99d0)then
C         mass = mass0
C      endif
C      if(mt.gt.99.99d0)then
C         mt = mt0
C      endif
*
      return
      end
***
