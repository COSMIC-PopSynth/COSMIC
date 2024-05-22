***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      integer kw,testflag
      real*8 lum,r,mt,mc,rl,z,teff,alpha
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap
      real*8 MLalpha
      external MLalpha
      parameter(lum0=7.0d+04,kap=-0.5d0)
*
*      windflag = 0 !BSE=0, startrack08=1, vink=2, vink+LBV for all
*      stars=3.
* Must be one of these values or mlwind will cause problem with code,
* i.e. mlwind not set (see last line of main if statement...).
    
      if(windflag.eq.0)then
* BSE
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
         dms = 0.d0
         if(lum.gt.4000.d0)then
            x = MIN(1.d0,(lum-4000.d0)/500.d0)
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
         endif
         if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            if(kw.gt.6)then
               dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
            else
               dms = MAX(dml,dms)
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dml,dms)
               endif
* LBV-like mass loss beyond the Humphreys-Davidson limit.
               x = 1.0d-5*r*sqrt(lum)
               if(lum.gt.6.0d+05.and.x.gt.1.d0)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
                  dms = dms + dml
               endif
            endif
         endif
*
         mlwind = dms
      elseif(windflag.eq.1)then
* StarTrack (Beclzynski+08)
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD, with no luminosity limit
* according to Belczynsk+08 pp. 174.
*
* This may not be what is actually assumed in StarTrack (see the windf1 function).
*
*
         dms = 0.d0
         if(lum.gt.4000.d0.or.(kw.ge.0.and.kw.le.1))then
            if(lum.gt.4000.d0)then
               x = MIN(1.d0,(lum-4000.d0)/500.d0)
            else
               x = 0.1d0/500.d0
            endif !or is it simply x = Min(1, lum/500)?
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
         endif
         if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            if(kw.gt.6)then
               dms = MAX(dml,1.0d-13*lum**(3.d0/2.d0)) !hewind here for KH06, not included for StarTrack...
            else
               dms = MAX(dml,dms)
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dml,dms)
               endif
* LBV-like mass loss beyond the Humphreys-Davidson limit.
               x = 1.0d-5*r*sqrt(lum)
               if(lum.gt.6.0d+05.and.x.gt.1.d0)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
                  dms = dms + dml
               endif
            endif
         endif
*
         mlwind = dms
      elseif(windflag.eq.2.or.windflag.eq.3)then
* Vink winds etc according to as implemented following
* Belczynski, Bulik, Fryer, Ruiter, Valsecchi, Vink & Hurley 2010.
*
* Firstly implement BSE 'old' winds that cover all other stars not
* accounted for by Vink winds (see Belczynski+09). Then implement
* Vink et al. winds.
*
* We also include the option for a variable metallicity-dependent mass
* loss parameter which eddlimflag is set, which makes the metallicity
* dependence become weaker as the star approaches the electron-scattering
* Eddington limit (Grafener & Hamann 2008, Giacobbo et al. 2018)
*
         teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
         dms = 0.d0
         if(lum.gt.4000.d0)then
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD after OB stars accounted for.
            x = MIN(1.d0,(lum-4000.d0)/500.d0)
            dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
            alpha = 0.5d0
            dms = dms*(z/zsun)**(alpha)
            testflag = 1
         endif
         if(kw.ge.2.and.kw.le.6)then
* 'Reimers' mass loss
            dml = neta*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
            if(kw.eq.5.or.kw.eq.6)then
               p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
               p0 = 10.d0**p0
               p0 = MIN(p0,2000.d0)
               dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
               dmt = 10.d0**dmt
               dmt = 1.d0*MIN(dmt,1.36d-09*lum)
               dml = MAX(dml,dmt)
            endif
            dms = MAX(dms,dml)
         endif
* Apply Vink, de Koter & Lamers (2001) OB star winds.
* Next check if hot massive H-rich O/B star in appropriate temperature ranges.
         if(teff.ge.12500.and.teff.le.25000)then
            if(eddlimflag.eq.0) alpha = 0.85d0
            if(eddlimflag.eq.1) alpha = MLalpha(mt,lum,kw)
            dms = -6.688d0 + 2.210d0*LOG10(lum/1.0d+05) -
     &            1.339d0*LOG10(mt/30.d0) - 1.601d0*LOG10(1.3d0/2.d0) +
     &            alpha*LOG10(z/zsun) + 1.07d0*LOG10(teff/2.0d+04)
            dms = 10.d0**dms
            testflag = 2
         elseif(teff.gt.25000.)then
*        Although Vink et al. formulae  are only defined until Teff=50000K,
*        we follow the Dutch prescription of MESA, and extend to higher Teff
             dms = -6.697d0 + 2.194d0*LOG10(lum/1.0d+05) -
     &            1.313d0*LOG10(mt/30.d0) - 1.226d0*LOG10(2.6d0/2.d0) +
     &            alpha*LOG10(z/zsun) +0.933d0*LOG10(teff/4.0d+04) -
     &            10.92d0*(LOG10(teff/4.0d+04)**2)
       dms = 10.d0**dms
       testflag = 2
         endif

         if((windflag.eq.3.or.kw.ge.2).and.kw.le.6)then
* LBV-like mass loss beyond the Humphreys-Davidson limit.
* Optional flag (windflag=3) to use for every non-degenerate star
* past the limit, rather than just for giant, evolved stars
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               if(eddlimflag.eq.0) alpha = 0.d0
               if(eddlimflag.eq.1) alpha = MLalpha(mt,lum,kw)
               dms = 1.5d0*1.0d-04*((z/zsun)**alpha)
               testflag = 3
            endif
         elseif(kw.ge.7.and.kw.le.9)then !WR (naked helium stars)
* If naked helium use Hamann & Koesterke (1998) WR winds reduced by factor of
* 10 (Yoon & Langer 2005), with Vink & de Koter (2005) metallicity dependence
            if(eddlimflag.eq.0) alpha = 0.86d0
            if(eddlimflag.eq.1) alpha = MLalpha(mt,lum,kw)
            dms = 1.0d-13*(lum**1.5d0)*((z/zsun)**alpha)
            testflag = 4
         endif
*
         mlwind = dms
      elseif(windflag.eq.4)then
*
* Calculate stellar wind mass loss following MIST as closely as possible
*
         dms = 0.d0
         if(kw.ge.0.and.kw.le.9)then
* 'Reimers' mass loss over the whole HRD
            dml = 0.1*4.0d-13*r*lum/mt
            if(rl.gt.0.d0) dml =
     &         dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Blocker 1995 for AGB.
            if(kw.eq.5.or.kw.eq.6)then
               dml = 4.83d-9*0.2*dml/0.1*((lum)**2.7)/((mt)**2.1)
            endif
            if(kw.gt.6)then
               dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
            else
               dms = MAX(dml,dms)
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dml,dms)
               endif
* LBV-like mass loss beyond the Humphreys-Davidson limit.
               x = 1.0d-5*r*sqrt(lum)
               if(lum.gt.6.0d+05.and.x.gt.1.d0)then
                  dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
                  dms = dms + dml
               endif
            endif
         endif
*
         mlwind = dms
      
      endif

      return
      end
***
