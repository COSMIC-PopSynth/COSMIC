***
      real*8 FUNCTION lzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Lzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      lzamsf = (a(1)*m**5*mx + a(2)*m**11)/
     &         (a(3) + m**3 + a(4)*m**5 + a(5)*m**7 +
     &          a(6)*m**8 + a(7)*m**9*mx)
*
      return
      end
***
      real*8 FUNCTION rzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Rzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      rzamsf = ((a(8)*m**2 + a(9)*m**6)*mx + a(10)*m**11 +
     &          (a(11) + a(12)*mx)*m**19)/
     &         (a(13) + a(14)*m**2 + 
     &          (a(15)*m**8 + m**18 + a(16)*m**19)*mx)
*
      return
      end
***
      real*8 FUNCTION tbgbf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the BGB or to
* Helium ignition if no FGB exists.
* (JH 24/11/97)
*
      tbgbf = (a(17) + a(18)*m**4 + a(19)*m**(11.d0/2.d0) + m**7)/
     &        (a(20)*m**2 + a(21)*m**7)
*
      return
      end
***
      real*8 FUNCTION tbgbdf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt mass.
* (JH 24/11/97)
*
      mx = SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*m**5*mx + m**7
      df = 4.d0*a(18)*m**3 + 5.5d0*a(19)*m**4*mx + 7.d0*m**6
      g = a(20)*m**2 + a(21)*m**7
      dg = 2.d0*a(20)*m + 7.d0*a(21)*m**6
      tbgbdf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION tbgdzf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt Z.
* (JH 14/12/98)
*
      mx = m**5*SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*mx + m**7
      df = a(117) + a(118)*m**4 + a(119)*mx
      g = a(20)*m**2 + a(21)*m**7
      dg = a(120)*m**2
      tbgdzf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION thookf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the end of the MS
* hook ( for those models that have one ) as a fraction of 
* the lifetime to the BGB
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      thookf = 1.d0 - 0.01d0*MAX(a(22)/m**a(23),a(24)+a(25)/m**a(26))
      thookf = MAX(thookf,0.5d0)
*
      return
      end
***
      real*8 FUNCTION ltmsf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the luminosity at the end of the MS
* (JH 24/11/97)
*
      ltmsf = (a(27)*m**3 + a(28)*m**4 + a(29)*m**(a(32)+1.8d0))/
     &        (a(30) + a(31)*m**5 + m**a(32))
* 
      return
      end
***
      real*8 FUNCTION lalphf(m)
      implicit none
      real*8 m,mcut,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity alpha coefficent.
* (JH 24/11/97)
*
      mcut = 2.d0
      if(m.ge.mcut)then
         lalphf = (a(33) + a(34)*m**a(36))/(m**0.4d0 + a(35)*m**1.9d0)
      else
         if(m.le.0.5d0)then
            lalphf = a(39)
         elseif(m.le.0.7d0)then
            lalphf = a(39) + ((0.3d0 - a(39))/0.2d0)*(m - 0.5d0)
         elseif(m.le.a(37))then
            lalphf = 0.3d0 + ((a(40)-0.3d0)/(a(37)-0.7d0))*(m - 0.7d0)
         elseif(m.le.a(38))then
            lalphf = a(40) + ((a(41)-a(40))/(a(38)-a(37)))*(m - a(37))
         else
            lalphf = a(41) + ((a(42)-a(41))/(mcut-a(38)))*(m - a(38))
         endif
      endif
*
      return
      end
***
      real*8 FUNCTION lbetaf(m)
      implicit none
      real*8 m,a1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity beta coefficent.
* (JH 24/11/97)
*
      lbetaf = a(43) - a(44)*m**a(45)
      lbetaf = MAX(lbetaf,0.d0)
      if(m.gt.a(46).and.lbetaf.gt.0.d0)then
         a1 = a(43) - a(44)*a(46)**a(45)
         lbetaf = a1 - 10.d0*a1*(m - a(46))
         lbetaf = MAX(lbetaf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION lnetaf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity neta exponent.
* (JH 24/11/97)
*
      if(m.le.1.d0)then
         lnetaf = 10.d0
      elseif(m.ge.1.1d0)then
         lnetaf = 20.d0
      else
         lnetaf = 10.d0 + 100.d0*(m - 1.d0)
      endif
      lnetaf = MIN(lnetaf,a(97))
*
      return
      end
***
      real*8 FUNCTION lhookf(m,mhook)
      implicit none
      real*8 m,mhook,a2,a(200)
      common /MSCFF/ a
*
* A function to evalute the luminosity at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         lhookf = 0.d0
      elseif(m.ge.a(51))then
         lhookf = MIN(a(47)/m**a(48),a(49)/m**a(50))
      else
         a2 = MIN(a(47)/a(51)**a(48),a(49)/a(51)**a(50))
         lhookf = a2*((m-mhook)/(a(51)-mhook))**0.4d0
      endif
*
      return
      end
***
      real*8 FUNCTION rtmsf(m)
      implicit none
      real*8 m,m2,rchk,a(200)
      common /MSCFF/ a
      real*8 rzamsf
      external rzamsf
*
* A function to evaluate the radius at the end of the MS
* Note that a safety check is added to ensure Rtms > Rzams
* when extrapolating the function to low masses. 
* (JH 24/11/97)
*
      m2 = a(62) + 0.1d0
      if(m.le.a(62))then
         rchk = 1.5d0*rzamsf(m)
         rtmsf = MAX(rchk,(a(52) + a(53)*m**a(55))/(a(54) + m**a(56)))
      elseif(m.ge.m2)then
         rtmsf = (a(57)*m**3+a(58)*m**a(61)+a(59)*m**(a(61)+1.5d0))/
     &           (a(60) + m**5)
      else
         rtmsf = a(63) + ((a(64) - a(63))/0.1d0)*(m - a(62))
      endif
* 
      return
      end
***
      real*8 FUNCTION ralphf(m)
      implicit none
      real*8 m,a5,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius alpha coefficent.
* (JH 24/11/97)
*
      if(m.le.0.5d0)then
         ralphf = a(73)
      elseif(m.le.0.65d0)then
         ralphf = a(73) + ((a(74) - a(73))/0.15d0)*(m - 0.5d0)
      elseif(m.le.a(70))then
         ralphf = a(74) + ((a(75)-a(74))/(a(70)-0.65d0))*(m - 0.65d0)
      elseif(m.le.a(71))then
         ralphf = a(75) + ((a(76) - a(75))/(a(71) - a(70)))*(m - a(70))
      elseif(m.le.a(72))then
         ralphf = (a(65)*m**a(67))/(a(66) + m**a(68))
      else
         a5 = (a(65)*a(72)**a(67))/(a(66) + a(72)**a(68))
         ralphf = a5 + a(69)*(m - a(72))
      endif
*
      return
      end
***
      real*8 FUNCTION rbetaf(m)
      implicit none
      real*8 m,m2,m3,b2,b3,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius beta coefficent.
* (JH 24/11/97)
*
      m2 = 2.d0
      m3 = 16.d0
      if(m.le.1.d0)then
         rbetaf = 1.06d0
      elseif(m.le.a(82))then
         rbetaf = 1.06d0 + ((a(81)-1.06d0)/(a(82)-1.d0))*(m-1.d0)
      elseif(m.le.m2)then
         b2 = (a(77)*m2**(7.d0/2.d0))/(a(78) + m2**a(79))
         rbetaf = a(81) + ((b2-a(81))/(m2-a(82)))*(m-a(82))
      elseif(m.le.m3)then
         rbetaf = (a(77)*m**(7.d0/2.d0))/(a(78) + m**a(79))
      else
         b3 = (a(77)*m3**(7.d0/2.d0))/(a(78) + m3**a(79))
         rbetaf = b3 + a(80)*(m - m3)
      endif
      rbetaf = rbetaf - 1.d0
*
      return
      end
***
      real*8 FUNCTION rgammf(m)
      implicit none
      real*8 m,m1,b1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius gamma coefficent.
* (JH 24/11/97)
*
      m1 = 1.d0
      if(m.gt.(a(88)+0.1d0))then
         rgammf = 0.d0
      else
         b1 = MAX(0.d0,a(83) + a(84)*(m1-a(85))**a(86))
         if(m.le.m1)then
            rgammf = a(83) + a(84)*ABS(m-a(85))**a(86)
         elseif(m.le.a(88))then
            rgammf = b1 + (a(89) - b1)*((m - m1)/(a(88) - m1))**a(87)
         else
            if(a(88).gt.m1) b1 = a(89)
            rgammf = b1 - 10.d0*b1*(m - a(88))
         endif
         rgammf = MAX(rgammf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION rhookf(m,mhook)
      implicit none
      real*8 m,mhook,m2,b2,a(200)
      common /MSCFF/ a
*
* A function to evalute the radius at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         rhookf = 0.d0
      elseif(m.le.a(94))then
         rhookf = a(95)*SQRT((m-mhook)/(a(94)-mhook))
      elseif(m.le.2.d0)then
         m2 = 2.d0
         b2 = (a(90) + a(91)*m2**(7.d0/2.d0))/
     &        (a(92)*m2**3 + m2**a(93)) - 1.d0
         rhookf = a(95) + (b2-a(95))*((m-a(94))/(m2-a(94)))**a(96)
      else
         rhookf = (a(90) + a(91)*m**(7.d0/2.d0))/
     &            (a(92)*m**3 + m**a(93)) - 1.d0
      endif
*
      return
      end
***
      real*8 FUNCTION lbgbf(m)
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the luminosity at the end of the 
* FGB ( for those models that have one )
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      lbgbf = (a(1)*m**a(5) + a(2)*m**a(8))/
     &        (a(3) + a(4)*m**a(7) + m**a(6))
* 
      return
      end
***
      real*8 FUNCTION lbgbdf(m)
      real*8 m,a(200)
      real*8 f,df,g,dg
      common /GBCFF/ a
*
* A function to evaluate the derivitive of the Lbgb function.
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      f = a(1)*m**a(5) + a(2)*m**a(8)
      df = a(5)*a(1)*m**(a(5)-1.d0) + a(8)*a(2)*m**(a(8)-1.d0)
      g = a(3) + a(4)*m**a(7) + m**a(6)
      dg = a(7)*a(4)*m**(a(7)-1.d0) + a(6)*m**(a(6)-1.d0)
*
      lbgbdf = (df*g - f*dg)/(g*g)
* 
      return
      end
***
      real*8 FUNCTION lbagbf(m,mhefl)
      implicit none
      real*8 m,mhefl,a4,a(200)
      common /GBCFF/ a
*
* A function to evaluate the BAGB luminosity. (OP 21/04/98)
* Continuity between LM and IM functions is ensured by setting
* gbp(16) = lbagbf(mhefl,0.0) with gbp(16) = 1.0.
*
      a4 = (a(9)*mhefl**a(10) - a(16))/(exp(mhefl*a(11))*a(16))
*
      if(m.lt.mhefl)then
         lbagbf = a(9)*m**a(10)/(1.d0 + a4*exp(m*a(11)))
      else
         lbagbf = (a(12) + a(13)*m**(a(15)+1.8d0))/(a(14) + m**a(15))
      endif
*
      return
      end
***
      real*8 FUNCTION rgbf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the GB.
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbf = a1*(lum**a(18) + a(17)*lum**a(19))
*
      return
      end
***
      real*8 FUNCTION rgbdf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the GB (as f(L)).
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &            a(17)*a(19)*lum**(a(19)-1.d0))
*
      return
      end
***
      real*8 FUNCTION ragbf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the AGB.
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbf = a1*(lum**a(18) + a(17)*lum**a4)
*
      return
      end
***
      real*8 FUNCTION ragbdf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the AGB (as f(L)).
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbdf = a1*(a(18)*lum**(a(18)-1.d0) + 
     &             a(17)*a4*lum**(a4-1.d0))
*
      return
      end
***
      real*8 FUNCTION mctmsf(m)
      implicit none
      real*8 m,m525
*
* A function to evaluate core mass at the end of the MS as a 
* fraction of the BGB value, i.e. this must be multiplied by 
* the BGB value (see below) to give the actual core mass (JH 5/9/99)
*
      m525 = m**(21.d0/4.d0)
      mctmsf = (1.586d0 + m525)/(2.434d0 + 1.02d0*m525)
*
      return
      end
***
      real*8 FUNCTION mcheif(m,mhefl,mchefl)
      implicit none
      real*8 m,mhefl,mchefl,mcbagb,a3,a(200)
      common /GBCFF/ a
      real*8 mcagbf
      external mcagbf
*
* A function to evaluate core mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars  (OP 25/11/97)
*
      mcbagb = mcagbf(m)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      mcheif = MIN(0.95d0*mcbagb,(a3 + a(33)*m**a(34))**(1.d0/4.d0))
*
      return
      end
***
      real*8 FUNCTION mheif(mc,mhefl,mchefl)
      implicit none
      real*8 mc,mhefl,mchefl,m1,m2,a3,a(200)
      common /GBCFF/ a
      real*8 mbagbf
      external mbagbf
*
* A function to evaluate mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars by inverting
* mcheif
*
      m1 = mbagbf(mc/0.95d0)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      m2 = ((mc**4 - a3)/a(33))**(1.d0/a(34))
      mheif = MAX(m1,m2)
*
      return
      end
***
      real*8 FUNCTION mcagbf(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate core mass at the BAGB  (OP 25/11/97)
*
      mcagbf = (a(37) + a(35)*m**a(36))**(1.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION mbagbf(mc)
      implicit none
      real*8 mc,mc4,a(200)
      common /GBCFF/ a
*
* A function to evaluate mass at the BAGB by inverting mcagbf.
*
      mc4 = mc**4
      if(mc4.gt.a(37))then
         mbagbf = ((mc4 - a(37))/a(35))**(1.d0/a(36))
      else
         mbagbf = 0.d0
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate Mc given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         mcgbtf = ((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (1.d0/(1.d0-GB(5)))
      else
         mcgbtf = ((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (1.d0/(1.d0-GB(6)))
      endif
*
      return
      end
***
      real*8 FUNCTION lgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate L given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         lgbtf = GB(4)*(((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (GB(5)/(1.d0-GB(5))))
      else
         lgbtf = GB(3)*(((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (GB(6)/(1.d0-GB(6))))
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbf(lum,GB,lx)
      implicit none
      real*8 lum,GB(10),lx
*
* A function to evaluate Mc given L for GB, AGB and NHe stars
*
      if(lum.le.lx)then
         mcgbf = (lum/GB(4))**(1.d0/GB(5))
      else
         mcgbf = (lum/GB(3))**(1.d0/GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lmcgbf(mc,GB)
      implicit none
      real*8 mc,GB(10)
*
* A function to evaluate L given Mc for GB, AGB and NHe stars
*
      if(mc.le.GB(7))then
         lmcgbf = GB(4)*(mc**GB(5))
      else
         lmcgbf = GB(3)*(mc**GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lHeIf(m,mhefl)
      implicit none
      real*8 m,mhefl,a(200)
      common /GBCFF/ a
*
* A function to evaluate He-ignition luminosity  (OP 24/11/97)
* Continuity between the LM and IM functions is ensured with a first
* call setting lhefl = lHeIf(mhefl,0.0)
*
      if(m.lt.mhefl)then
         lHeIf = a(38)*m**a(39)/(1.d0 + a(41)*EXP(m*a(40)))
      else
         lHeIf = (a(42) + a(43)*m**3.8d0)/(a(44) + m**2)
      endif
*
      return
      end
***
      real*8 FUNCTION lHef(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the ratio LHe,min/LHeI  (OP 20/11/97)
* Note that this function is everywhere <= 1, and is only valid
* for IM stars
*
      lHef = (a(45) + a(46)*m**(a(48)+0.1d0))/(a(47) + m**a(48))
*
      return
      end
***
      real*8 FUNCTION rminf(m)
      implicit none
      real*8 m,mx,a(200)
      common /GBCFF/ a
*
* A function to evaluate the minimum radius during He-burning
* for IM & HM stars  (OP 20/11/97)
*
      mx = m**a(53)
      rminf = (a(49)*m + (a(50)*m)**a(52)*mx)/(a(51) + mx)
*
      return
      end
***
      real*8 FUNCTION tHef(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a(200)
      common /GBCFF/ a
      real*8 themsf
      external themsf
*
* A function to evaluate the He-burning lifetime.  (OP 26/11/97)
* For IM & HM stars, tHef is relative to tBGB.
* Continuity between LM and IM stars is ensured by setting
* thefl = tHef(mhefl,0.0,,0.0), and the call to themsf ensures
* continuity between HB and NHe stars as Menv -> 0.
*
      if(m.le.mhefl)then
         mm = MAX((mhefl - m)/(mhefl - mc),1.0d-12)
         tHef = (a(54) + (themsf(mc) - a(54))*mm**a(55))*
     &          (1.d0 + a(57)*EXP(m*a(56)))
      else
         mm = m**5
         tHef = (a(58)*m**a(61) + a(59)*mm)/(a(60) + mm)
      endif
*
      return
      end
***
      real*8 FUNCTION tblf(m,mhefl,mfgb)
      implicit none
      real*8 m,mhefl,mfgb,mr,m1,m2,r1,a(200)
      common /GBCFF/ a
      real*8 lheif,rminf,ragbf
      external lheif,rminf,ragbf
*
* A function to evaluate the blue-loop fraction of the He-burning
* lifetime for IM & HM stars  (OP 28/01/98)
*
      mr = mhefl/mfgb
      if(m.le.mfgb) then
         m1 = m/mfgb
         m2 = log10(m1)/log10(mr)
         m2 = max(m2,1.0d-12)
         tblf = a(64)*m1**a(63) + a(65)*m2**a(62)
      else
         r1 = 1.d0 - rminf(m)/ragbf(m,lheif(m,mhefl),mhefl)
         r1 = max(r1,1.0d-12)
         tblf = a(66)*m**a(67)*r1**a(68)
      end if
      tblf = MIN(1.d0,MAX(0.d0,tblf))
      if(tblf.lt.1.0d-10) tblf = 0.d0
*
      return
      end
***
      real*8 FUNCTION lzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a4,a5,a(200)
      common /GBCFF/ a
      real*8 lzhef
      external lzhef
*
* A function to evaluate the ZAHB luminosity for LM stars. (OP 28/01/98)
* Continuity with LHe,min for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to lzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      a5 = lzhef(mc)
      a4 = (a(69) + a5 - a(74))/((a(74) - a5)*exp(a(71)*mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      lzahbf = a5 + (1.d0 + a(72))*a(69)*mm**a(70)/
     &         ((1.d0 + a(72)*mm**a(73))*(1.d0 + a4*EXP(m*a(71))))
*
      return
      end
***
      real*8 FUNCTION rzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,rx,ry,mm,f,a(200)
      common /GBCFF/ a
      real*8 rzhef,rgbf,lzahbf
*
* A function to evaluate the ZAHB radius for LM stars. (OP 28/01/98)
* Continuity with R(LHe,min) for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to rzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      rx = rzhef(mc)
      ry = rgbf(m,lzahbf(m,mc,mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      f = (1.d0 + a(76))*mm**a(75)/(1.d0 + a(76)*mm**a(77))
      rzahbf = (1.d0 - f)*rx + f*ry
*
      return
      end
***
      real*8 FUNCTION lzhef(m)
      implicit none
      real*8 m,m15
*
* A function to evaluate Naked Helium star 'ZAMS' luminosity
*
      m15 = m*SQRT(m)
      lzhef = 1.5262d+04*m**(41.d0/4.d0)/
     &        (0.0469d0 + m**6*(31.18d0 + m15*(29.54d0 + m15)))
*
      return
      end
***
      real*8 FUNCTION rzhef(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star 'ZAMS' radius
*
      rzhef = 0.2391d0*m**4.6d0/(0.0065d0 + (0.162d0 + m)*m**3)
*
      return
      end
***
      real*8 FUNCTION themsf(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star main sequence lifetime
*
      themsf = (0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6)/m**(13.d0/2.d0)
*
      return
      end
***
      real*8 FUNCTION rhehgf(m,lum,rx,lx)
      implicit none
      real*8 m,lum,rx,lx,cm
*
* A function to evaluate Helium star radius on the Hertzsprung gap 
* from its mass and luminosity. 
*
      cm = 2.0d-03*m**(5.d0/2.d0)/(2.d0 + m**5)
      rhehgf = rx*(lum/lx)**0.2d0 + 0.02d0*(EXP(cm*lum) - EXP(cm*lx))
*
      return
      end
***
      real*8 FUNCTION rhegbf(lum)
      implicit none
      real*8 lum
*
* A function to evaluate Helium star radius on the giant branch. 
*
      rhegbf = 0.08d0*lum**(3.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION lpertf(m,mew)
      implicit none
      real*8 m,mew
      real*8 b,c
*
* A function to obtain the exponent that perturbs luminosity.
*
      b = 0.002d0*MAX(1.d0,2.5d0/m)
      c = 3.d0
      lpertf = ((1.d0 + b**c)*((mew/b)**c))/(1.d0+(mew/b)**c)
*
      return
      end
***
      real*8 FUNCTION rpertf(m,mew,r,rc)
      implicit none
      real*8 m,mew,r,rc
      real*8 a,b,c,q,fac,facmax
*
* A function to obtain the exponent that perturbs radius.
*
      if(mew.le.0.d0)then
         rpertf = 0.d0
      else
         a = 0.1d0
         b = 0.006d0*MAX(1.d0,2.5d0/m)
         c = 3.d0
         q = log(r/rc)
         fac = a/q
         facmax = -14.d0/log10(mew)
         fac = MIN(fac,facmax)
         rpertf = ((1.d0 + b**c)*((mew/b)**c)*(mew**fac))/
     &            (1.d0+(mew/b)**c)
      endif
*
      return
      end
***
      real*8 FUNCTION vrotf(m,ST)
      implicit none
      integer ST
      real*8 m
*
      if(ST.gt.0)then
         if(m.gt.6.35d0)then
            vrotf = (10.d0*m**(-0.0354d0))/(0.0389d0+m**(-7.95d0))
         else
            vrotf = (13.4d0*m**(-0.12d0))/(0.0389d0+m**(-7.95d0))
         endif
      else
         vrotf = 330.d0*m**3.3d0/(15.d0 + m**3.45d0)
      endif
*
      return
      end
***
      real*8 FUNCTION celamf(kw,m,lum,rad,rzams,menv,fac)
      implicit none
      integer kw
      real*8 m,lum,rad,rzams,menv,fac
      real*8 lam1,lam2,m1,logm,logl
      real*8 aa,bb,cc,dd
*
* A function to estimate lambda for common-envelope.
*
      if(fac.ge.0.d0)then
*
* No fits yet for naked He stars...
*
         if(kw.gt.6)then
            celamf = 0.5d0
            goto 90
         endif
*
         if(menv.gt.0.d0)then
* Formulae for giant-like stars; also used for HG and CHeB stars close
* to the Hayashi track.
            logl = log10(lum)
            logm = log10(m)
            if(kw.le.5)then
               m1 = m
               if(kw.gt.3) m1 = 100.d0
               lam1 = 3.d0/(2.4d0 + 1.d0/m1**(3.d0/2.d0)) - 0.15d0*logl
               lam1 = MIN(lam1,0.8d0)
            else
               lam1 = -3.5d0 - 0.75d0*logm + logl
            endif
            if(kw.gt.3)then
               lam2 = MIN(0.9d0,0.58d0 + 0.75d0*logm) - 0.08d0*logl
               if(kw.lt.6)then
                  lam1 = MIN(lam2,lam1)
               else
                  lam1 = MAX(lam2,lam1)
                  lam1 = MIN(lam1,1.d0)
               endif
            endif
            lam1 = 2.d0*lam1
            if(fac.gt.0.d0)then
* Use a fraction FAC of the ionization energy in the energy balance.
               if(kw.le.3)then
                  aa = MIN(1.2d0*(logm - 0.25d0)**2 - 0.7d0,-0.5d0)
               else
                  aa = MAX(-0.2d0 - logm,-0.5d0)
               endif
               bb = MAX(3.d0 - 5.d0*logm,1.5d0)
               cc = MAX(3.7d0 + 1.6d0*logm,3.3d0 + 2.1d0*logm)
               lam2 = aa + ATAN(bb*(cc - logl))
               if(kw.le.3)then
                  dd = MAX(0.d0,MIN(0.15d0,0.15d0 - 0.25d0*logm))
                  lam2 = lam2 + dd*(logl - 2.d0)
               endif
               lam2 = MAX(lam2,1.d-2)
               lam2 = MAX(1.d0/lam2,lam1)
               if(fac.ge.1.d0)then
                  lam1 = lam2
               else
                  lam1 = lam1 + fac*(lam2 - lam1)
               endif
            endif
         endif
*
         if(menv.lt.1.d0)then
* Formula for HG stars; also reasonable for CHeB stars in blue loop.
            lam2 = 0.42d0*(rzams/rad)**0.4d0
* Alternatively:
*           lam2 = 0.3d0*(rtms/rad)**0.4d0
            lam2 = 2.d0*lam2
         endif
*
         if(menv.le.0.d0)then
            celamf = lam2
         elseif(menv.ge.1.d0)then
            celamf = lam1
         else
* Interpolate between HG and GB values depending on conv. envelope mass.
            celamf = lam2 + sqrt(menv)*(lam1 - lam2)
         endif
*
 90      continue
*
      else
         celamf = -1.d0*fac
      endif
*
      return
      end
***

      real*8 FUNCTION celamf_xu_li(kw,m0,m,mc,R,fac)
      implicit none
      integer kw,GBt,nRGBbc,ig,ii
      real*8 m0,m,mc,lum,R,rzams,menv,fac
      real*8 logM,logR,LMHMbs,calc_logBE,lamda
      real*8 logBE0,Mzams,logRbd,dlogBE
      integer ndat(6),ms(325),rs(325)
      double precision alphas(325) 
      real*8 RGBb_coef(5)
      PARAMETER(logBE0 = 33.29886d0)
      ndat = (/0,8,44,144,309,325/)

      logM = log10(M)
      logR = log10(R)
      nRGBbc = 5
      Mzams = m0
      ms = (/ 0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,
     &2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,
     &5,5,5,5,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,
     &1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
     &3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,
     &5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,
     &7,7,7,7,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,
     &9,9,9,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
     &1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,
     &2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,
     &3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,
     &5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,
     &6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,
     &7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,
     &9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,
     &10,10,10,10,10,10,10,10,10,0,0,0,0,1,1,1,1,2,2,2,
     &2,3,3,3,3 /)
      rs = (/0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5,
     &0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,
     &2,3,4,5,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,
     &6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,
     &6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,
     &6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,
     &6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,
     &6,7,8,9,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,
     &-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,-3,-2,-1,0,1,
     &2,3,4,5,6,7,8,9,10,-4,-3,-2,-1,0,1,2,3,4,5,6,
     &7,8,9,10,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,
     &-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,-3,-2,-1,0,1,
     &2,3,4,5,6,7,8,9,10,-4,-3,-2,-1,0,1,2,3,4,5,6,
     &7,8,9,10,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,
     &-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,-4,-3,-2,-1,0,1,
     &2,3,4,5,6,7,8,9,10,0,1,2,3,0,1,2,3,0,1,2,
     &3,0,1,2,3/)
      alphas = (/1.50305864307172658556d+01,
     &4.96110465703095582235d-01,
     &-9.16608282566077403608d-01,2.48033316991016244968d-01,
     &1.77238406649739155263d+00,-6.52669943496246296455d-01,
     &6.23616025599144307989d-01,-1.77703681033998805994d-01,
     &1.56617552673034907684d+01,-3.38271409511212528543d+00,
     &3.99712000854973803499d+00,-3.02506056136319223526d+00,
     &1.12651646333780397491d+00,-1.62914377442334834534d-01,
     &7.36910246144982838956d+00,-2.30880577686855161801d+01,
     &3.32294193400783512971d+01,-2.22327077943441828722d+01,
     &6.85986578867425578210d+00,-7.85596865167521363205d-01,
     &3.04288311387925105578d+01,-1.38697534163042433875d+02,
     &2.72764619820308155340d+02,-2.82250602598139892052d+02,
     &1.37077309352478550863d+02,-2.43363801583426564434d+01,
     &-2.40391883810966504598d+02,7.20123937632566253342d+02,
     &-9.56047390568349555906d+02,8.80126339460319968566d+02,
     &-4.43060386994394036719d+02,8.34534841118213392974d+01,
     &1.19040215556355087756d+03,-3.06320152050599153881d+03,
     &2.78453473219604893529d+03,-1.43582241909097660937d+03,
     &4.99720580375658244066d+02,-8.44179485778146130315d+01,
     &-1.27279708339606577283d+03,3.68563309227265335721d+03,
     &-3.66200350857788362191d+03,1.74828529177925724980d+03,
     &-4.43446115665101842751d+02,5.16038082522078980219d+01,
     &-1.51049254914130983707d+03,5.90804592923560267081d+03,
     &-9.78411552163013584504d+03,9.03281646038073631644d+03,
     &-5.04733341717106668511d+03,1.70840439977420237483d+03,
     &-3.16251664756395882705d+02,1.74202707579704139107d+01,
     &3.67713293499369164863d+00,-5.09055045574688058707d-01,
     &1.18642296961242536781d+05,-4.22882345998080738354d+05,
     &6.02728568485423107632d+05,-4.09618748060796875507d+05,
     &8.90979215331482701004d+04,5.95414095173196401447d+04,
     &-5.18590927506359803374d+04,1.72331732723822569824d+04,
     &-2.81277261133411593619d+03,1.86487606226860776815d+02,
     &-1.40854420295521779917d+06,3.16920022432578727603d+06,
     &4.61300052596119523514d+05,-8.61105917229603044689d+06,
     &1.23556398056155946106d+07,-8.94036121833491511643d+06,
     &3.81112897666227677837d+06,-9.71158692551792715676d+05,
     &1.37401910653425060445d+05,-8.32876128728674302693d+03,
     &1.72458181532576158643d+07,-4.52772539129179418087d+07,
     &2.01148809934781342745d+07,6.09227038679622411728d+07,
     &-1.07679655401528432965d+08,8.29801717072715312243d+07,
     &-3.64493890278426855803d+07,9.44902866139182634652d+06,
     &-1.35153547941456316039d+06,8.25342644314703647979d+04,
     &-1.36692876826326847076d+08,4.51025336598177015781d+08,
     &-5.31011502233843207359d+08,1.64806847028889596462d+08,
     &2.15493288928440541029d+08,-2.68977365806154668331d+08,
     &1.38776598747677832842d+08,-3.89913542043956518173d+07,
     &5.84939274544744472951d+06,-3.68464314099946233910d+05,
     &5.88007632627574801445d+08,-2.19546849066643857956d+09,
     &3.29578755621284484863d+09,-2.44640117334454774857d+09,
     &7.63562415504432678223d+08,1.37457055430822163820d+08,
     &-2.06529556624937295914d+08,7.50013034409245699644d+07,
     &-1.26735707259846013039d+07,8.54653687692407402210d+05,
     &-1.38327254361235642433d+09,5.60614067713224506378d+09,
     &-9.43342679709010314941d+09,8.54102479910387897491d+09,
     &-4.42476607837071800232d+09,1.22015302616584563255d+09,
     &-9.17636537403069436550d+07,-4.19402333332933112979d+07,
     &1.20023414749641995877d+07,-9.83676207032694481313d+05,
     &1.70433842457238912582d+09,-7.46564615829306602478d+09,
     &1.36213117744959812164d+10,-1.36449069424619369507d+10,
     &8.23382455328104782104d+09,-3.05178841976632213593d+09,
     &6.64362036529520750046d+08,-7.09507497840798497200d+07,
     &7.89097510117984027602d+05,3.47019468185705947690d+05,
     &-9.03133463246678471565d+08,4.47621684993862724304d+09,
     &-8.96349174100218391418d+09,9.77159178772765159607d+09,
     &-6.46222307786082077026d+09,2.69448710911580371857d+09,
     &-7.04153944276203036308d+08,1.09216364710150659084d+08,
     &-8.71251345965249091387d+06,2.35094345904298679670d+05,
     &7.47826861721657812595d+07,-6.76538905537033677101d+08,
     &1.73090915363686656952d+09,-2.18473390758292007446d+09,
     &1.61084727759012627602d+09,-7.39637661531803250313d+08,
     &2.13872151698962450027d+08,-3.76329373731744587421d+07,
     &3.63438629714747751132d+06,-1.44132561512217595009d+05,
     &2.53729226092276048660d+09,-9.90007540851090812683d+09,
     &1.52593834113249797821d+10,-1.01549528384376564026d+10,
     &2.99829578010063916445d+07,5.78893490124735832214d+09,
     &-6.1460989849077844619d+09,4.04278241615104389191d+09,
     &-1.80654362590145325661d+09,4.91334348953023314476d+08,
     &-5.44061721887778788805d+07,-6.47383359300961066037d+06,
     &1.76558650155847985297d+06,1.00096546886151510989d+05,
     &-3.39961914804177431506d+04,-1.43765508369811935425d+10,
     &5.53386381611402206421d+10,-8.64972441405090179443d+10,
     &6.18263781388091888428d+10,-1.27031002818969459534d+10,
     &-1.22741674009294738770d+10,1.32545543793032627106d+10,
     &-7.72103742705006504059d+09,2.99797494889310359955d+09,
     &-4.51316994491591095924d+08,-2.12359218253416419029d+08,
     &1.28106564253871873021d+08,-2.52927841362378038466d+07,
     &1.52832625190069107339d+06,5.33991439763792805024d+04,
     &3.53888765741276016235d+10,-1.30949004957460281372d+11,
     &2.03740123625964263916d+11,-1.45485633356229827881d+11,
     &3.38169397075351562500d+10,1.75117599864639778137d+10,
     &-1.68526525631832427979d+10,7.96671769498775196075d+09,
     &-3.11421864218404531479d+09,6.72407024087293863297d+08,
     &1.69479333713295638561d+08,-1.48711963495587915182d+08,
     &3.25634279867792949080d+07,-1.75731226257729995996d+06,
     &-1.35079742770175478654d+05,-5.02259172277233886719d+10,
     &1.70340081645104858398d+11,-2.59178892024369842529d+11,
     &1.80456668127881774902d+11,-4.10109717020726623535d+10,
     &-1.59603847360727577209d+10,1.16582282107490043640d+10,
     &-2.69128104888696050644d+09,7.08683177906285405159d+08,
     &-4.34911479088057935238d+08,1.01044743218226939440d+08,
     &1.50179933915116582066d+07,-4.86460542015785723925d+06,
     &-1.18190563420071476139d+06,2.78623554417385661509d+05,
     &4.64705908425810546875d+10,-1.30875109087104644775d+11,
     &1.84384613263450225830d+11,-1.14315528459969711304d+11,
     &1.54155640604285240173d+10,1.32641401056199913025d+10,
     &-5.51530197341053485870d+09,-1.70800161550766736269d+08,
     &7.78830644995285868645d+08,-2.06842292132320821285d+08,
     &-2.58832894926396012306d+07,2.22092620407946482301d+07,
     &-6.95661957240117993206d+06,1.92690144883910729550d+06,
     &-2.33277229786988784326d+05,-3.03660798334242973328d+10,
     &5.92434377201234664917d+10,-6.54677950676404571533d+10,
     &1.98006473040911102295d+10,1.80741102193409042358d+10,
     &-1.24110750642199649811d+10,2.06936088729867076874d+09,
     &-5.55200935775747060776d+08,3.31793995932478666306d+08,
     &-3.77282162409072965384d+07,3.77271593664872169029d+05,
     &-2.39515931142189726233d+06,4.78976734628147096373d+05,
     &-2.24492564827198133571d+05,4.79792411042873136466d+04,
     &1.39020977034095058441d+10,-1.24214230056445674896d+10,
     &3.40159755029327201843d+09,1.62527040079761524200d+10,
     &-2.40572617307258148193d+10,1.01576765639793376923d+10,
     &-2.79329356560366511345d+08,-3.70225353787005186081d+08,
     &-1.05583829363182082772d+08,4.22791698590517193079d+07,
     &-2.37510881991530815139d+06,-2.48551971484777750447d+06,
     &2.30980219498265953735d+06,-5.97702223240277264267d+05,
     &4.54666582402248095605d+04,-2.83252506215792465210d+09,
     &-5.8710759361964426040d+09,9.35640963237767410278d+09,
     &-9.10900080282122421265d+09,9.58587735224136924744d+09,
     &-5.49030619645441913605d+09,1.36275859178348875046d+09,
     &-8.57005601597409099340d+07,-7.86904883546570241451d+07,
     &4.99692720418420732021d+07,-7.09569033777875266969d+06,
     &-8.55338083130266284570d+05,-5.61264059584917849861d+05,
     &3.08073506627650989685d+05,-3.28288523737060459098d+04,
     &-1.46525129680925726891d+09,9.31328032341473007202d+09,
     &-1.05576564874556827545d+10,4.15304574768818807602d+09,
     &-5.00505353118861734867d+08,2.38979281271329969168d+08,
     &-3.43075090044980287552d+08,1.53399431073772042990d+08,
     &6.15683430398282781243d+06,-2.10229696774827726185d+07,
     &6.68991617223582346924d+05,2.39900342987180408090d+06,
     &-4.43480899599777883850d+05,-2.42983373542564295349d+04,
     &7.76394675149368777056d+03,1.21311579654507946968d+09,
     &-5.24070703990287113190d+09,7.06700234139572620392d+09,
     &-4.37710070889741706848d+09,1.25442565132006978989d+09,
     &-4.98846050957186818123d+07,-6.82854203030566722155d+07,
     &3.25620574902522787452d+07,-2.20218804932029694319d+07,
     &9.02162431089685671031d+06,-2.40492887758824275807d+05,
     &-8.68585806973894243129d+05,2.25402435254328796873d+05,
     &-1.55503066541077005240d+04,-3.94672900362825259890d+02,
     &-2.52924874342063546181d+08,1.08433769070899581909d+09,
     &-1.75274841405356764793d+09,1.53857146715462470055d+09,
     &-8.55304471302122354507d+08,3.31495823209429860115d+08,
     &-9.56754638305779546499d+07,1.92680795062012895942d+07,
     &-1.01626479920090991072d+06,-6.69255818740246933885d+05,
     &1.56196811545844484499d+04,9.97457630392971332185d+04,
     &-2.92981386761484645831d+04,2.94914409981776043423d+03,
     &-6.54138557526380992613d+01,1.34142729153601880654d+01,
     &-8.33977084929894418863d-01,6.31803234438807592710d-01,
     &-1.70217751115103232973d-01,1.57297096125805579980d+00,
     &3.04879789211891349954d-01,-6.10228805816398045536d-01,
     &2.31090540934666771600d-01,-1.41108730414326188907d+00,
     &1.23004856793249794933d+00,-2.83350980192163703908d-01,
     &-4.69140139332027000796d-02,5.62788156340448986192d-01,
     &-6.89182559641456582433d-01,2.69502539954297459790d-01,
     &-2.14046891159626988255d-02/) 
      RGBb_coef = (/2.406370d-01,1.089220d+00,1.953180d+00,
     &-1.032180d+00,0.000000d+00/)



*
      if(fac .ge. 0.d0) then
       
         if((kw .gt. 6) .or. (kw.le.1) ) then
            celamf_xu_li = 0.5d0
            goto 90
         endif
      endif
*
      if((kw .eq. 2) .or. (kw .eq. 3)) then
          GBt = 1
      else
          GBt = 2
      endif
*
      LMHMbs = log10(11.7d0)
   
*     Find the group from the mass and evolutionary state (ig = 1-4):
      if(logM.le.LMHMbs) then                           
*        group LM (low mass)
         if(GBt.eq.1) then
*           group LMR (low-mass RGB)
            ig = 1

*           Compute boundary radius between LMR1 an LMR2:
            logRbd = 0.0d0
            do 10 ii=1,nRGBbc
               logRbd=logRbd+RGBb_coef(ii)*logM**(ii-1)
10          continue
            if(logR.gt.logRbd) then
*              group LMR2 (low-mass RGB 2)
               ig = 2 
         else
*           group LMA (low-mass AGB)
            ig = 3
         end if
      else
*        group HM (high mass)
         ig = 4
      end if
*
*     Compute the 10-logarithm of the binding energy:
      calc_logBE = 0.0d0
*     Avoid 0**0
      do 11 ii=ndat(ig)+1,ndat(ig+1)
         dlogBE = alphas(ii)*(logM+tiny(logM))**dble(ms(ii))* 
     &  (logR+tiny(logR))**dble(rs(ii))
         calc_logBE = calc_logBE + dlogBE

11    continue
      
      lamda = 1.0d0 + 0.25d0 * ((Mzams-M)/Mzams)**2
*
      calc_logBE = lamda*calc_logBE
   
   
*     BE was originally expressed in erg/solar mass to avoid large numbers, so we need to convert to erg here:
      calc_logBE = calc_logBE + logBE0
*
      celamf_xu_li = 3.794d41*m*(m-mc)/R/10**(calc_logBE-7.0d0)
*
      if(abs(celamf_xu_li) .lt. 0.001d0) celamf_xu_li = 0.001d0
90       continue
      else
         celamf_xu_li = -1.0*fac
      endif
      return
      END

