***
      SUBROUTINE hrdiag_remnant(zpars,mt,mc,lum,r,aj,kw)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 zpars(20)

      real*8 mt,mc,lum,r,aj
      integer kw

      real*8 mch,mchold
      real*8 xx,fac,ahe,aco
      parameter(ahe=4.d0,aco=16.d0)

      
* input mc(or mcmax),mass
* output mt, kw, mc, mass,bhspin
* common: mch, mxns,
* ifflag,ecsn (or mc1), ecsn_low(mc2), bhspinflag,
*  bhspinmag,  Mbh_initial

         mch = 1.44d0 !set here owing to AIC ECSN model.

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
*                  mt = 1.38d0
                  if(wd_mass_lim.eq.1)then
                      mt = 1.38d0
                  else
                      mt = mc
                  endif
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
* Store the initial BH mass for calculating the ISCO later
         if(Mbh_initial.eq.0)then
            Mbh_initial = mt
         endif
         lum = 1.0d-10
         r = 4.24d-06*mt
      endif
*
      !added by PA
      if(kw .eq.15) then
*
*         Massless remnant
*
          aj = 0.d0
          mt = 0.d0
          lum = 1.0d-10
          r = 1.0d-10
      endif
*
      end
