***
* New flag introduced (August 2023)
* rtmsflag = 0 uses the SSE rtms for stars EXCEPT if M>200Msun 
* and Z<0.0008(~0.04Zsun). For the exception, an extrapolation
* is used (check the code in function rtmsf for details). This 
* ad hoc extrapolation works well for stars with M < 4000Msun 
* for Z >= 0.01Zsun. For lower metallicities, you may still run
* into issues of rtms<0 for very massive stars (i.e. for Z<0.01Zsun
* use at your own risk)
*
* rtmsflag > 0 uses rtms from simulation data (1=[BoOST], 2=[BPASS])
* Only common metallicity cases between BoOST & BPASS are coded.
*
* BoOST metallicities = [1.1e-4, 2.1e-4, 1e-3, 2e-3] corresponding
* to [dwarfD, IZw18, dwarfA, SMC] models (Szecsi et al. (2022))
*
* BPASS metallicities = [1e-4, *2e-4*, 1e-3, 2e-3]. 
* NOTE : For BPASS, we used a power law to fit the rtms v/s mzams 
* values for each metallicity. We have coded the best fit power-laws
* for the above BPASS metallicities here.Since, Z = 2e-4 model is 
* not available in the BPASS tracks, we assume the same rtms v/s 
* mzams power law as Z=1e-4 for Z=2e-4.
*
***
      real*8 FUNCTION lzamsf(m)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
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
      real*8 FUNCTION rtmssse(m)
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
         rtmssse = MAX(rchk,(a(52) + a(53)*m**a(55))/(a(54) + m**a(56)))
      elseif(m.ge.m2)then
         rtmssse = (a(57)*m**3+a(58)*m**a(61)+a(59)*m**(a(61)+1.5d0))/
     &           (a(60) + m**5)
      else
         rtmssse = a(63) + ((a(64) - a(63))/0.1d0)*(m - a(62))
      endif
*
      return
      end

***
*** NEW CODE BETWEEN THESE LINES (BEGIN)

      real*8 FUNCTION rtmsf(m, met)
      implicit none
      INCLUDE 'const_bse.h'
      real*8 m,met,Rtms200,Rtms199,slope
      real*8 rtmssse, rtmsBoost, rtmsBpass
      external rtmssse
*
* For Z < 0.04 Z_sun (i.e. Z < 0.0008) the polynomial fitting extrapolation of
* Hurley et al. leads to negative stellar radii at the end of the main sequence.
* This fix calculates the slope of the Rtms-M curve for a given Z at M=200 and
* extrapolates the Rtms linearly for M>200Msun assuming the same slope. The fix
* is only applied to stars with Z < 0.0008 and M > 200Msun
*
      if (rtmsflag .eq. 0) then
*         print*,'inside cosmic conditional'
         if( (m .le. 200.0) .OR. (met .ge. 0.0008) )then
            rtmsf = rtmssse(m)
         else
            Rtms200 = rtmssse(DFLOAT(200))
            Rtms199 =  rtmssse(DFLOAT(199))
            slope = Rtms200 - Rtms199;
            rtmsf = Rtms200 + slope*(m-200.0)
*            print*,'RTMSF inside the else', met,m,Rtms200, Rtms199,
*     &               slope,rtmsf
         endif

      elseif (rtmsflag .eq. 1) then
*         print*,'inside Boost conditional'
         rtmsf = rtmsBoost(m, met)

      elseif (rtmsflag .eq. 2) then
*         print*,'inside BPASS conditional'
         rtmsf = rtmsBpass(m, met)

      else
         print*,'Specify a valid value of rtmsflag'

      endif
*      print*,'met, m, rtmsf, rtmsflag', met, m, rtmsf, rtmsflag
      return
      end
*
      real*8 FUNCTION interp(val, arrsize, arr1, arr2, extrapolate)
          implicit none
          character(*) extrapolate
          integer j, idx_u, idx_l, arrsize
          real*8 val, arr1(arrsize), arr2(arrsize), num, den, slp

          if (arr1(1) .gt. val) then
              idx_u = 2
              idx_l = 1
              num = (arr2(idx_u) - arr2(idx_l))
              den = arr1(idx_u) - arr1(idx_l)
              slp = num/den
              interp = arr2(idx_l) - slp * (arr1(idx_l) - val)
*              print*,'inside 1st if', interp, arr1(1), val
          else
              if (arr1(arrsize) .ge. val) then
                  do j = 1, arrsize
                      if (arr1(j) .ge. val) then
*                      print*, 'arr1[j]', arr1(j)
                      idx_u = j
                      exit
                      end if
                  end do
                  idx_l = idx_u - 1
                  num = (arr2(idx_u) - arr2(idx_l))
                  den = arr1(idx_u) - arr1(idx_l)
                  slp = num/den
                  interp = slp * (val - arr1(idx_l)) + arr2(idx_l)
*                  print*,'inside else - if', num, den
              else
                  if (extrapolate == 'const') then
                      interp = arr2(arrsize)
*                      print*,'inside else - else - if'
                  else
                      idx_u = arrsize
                      idx_l = idx_u - 1
                      num = (arr2(idx_u) - arr2(idx_l))
                      den = arr1(idx_u) - arr1(idx_l)
                      slp = num/den
                      interp = slp*(val - arr1(idx_u)) + arr2(idx_u)
*                      print*,'inside else - else - else'
                  endif
              endif
          endif
      return
      end FUNCTION interp
****
      real*8 FUNCTION rtmsBoost(m, met)
      implicit none
      external interp
      real*8 m,met,met_cut1,met_cut2
      real*8 interp, rzams, rtms
      real*8 mBoost(10), RZBoost(10), RTBoost(10)
      integer j, idx_l, idx_u
      met_cut1 = 0.9*met
      met_cut2 = 1.1*met
*      print*, 'mass is ', m, LOG10(m)
*      print*,'met_cut1, met_cut2 :', met_cut1, met_cut2
*
      if ( (1.1d-4 .ge. met_cut1) .and. (1.1d-4 .le. met_cut2) ) then
         mBoost = (/9., 12., 19., 30., 40., 55., 80., 100.,250.,560./)
         RZBoost = (/2.43392, 2.92073, 3.86161, 4.99635, 5.86427,
     &      6.98092, 8.5704 , 9.68978, 16.1026 , 24.6849/)
         RTBoost = (/6.68652, 8.35836, 12.2768, 19.1216, 26.1014,
     &      34.7308, 50.2973, 63.3516, 696.887, 477.342/)
*         print*,'inside if'
*
      elseif ((2.1d-4 .ge. met_cut1) .and. (2.1d-4 .le. met_cut2)) then
         mBoost = (/9.0d0, 1.7d1, 2.0d1, 2.6d1, 4.5d1,
     &     7.7d1, 1.0d2, 1.5d2, 2.57d2, 5.75d2/)
         RZBoost = (/2.55385, 3.78012, 4.16676, 4.84233, 6.56514,
     &      8.81132, 10.1899, 12.7489, 17.242, 26.9047/)
         RTBoost = (/7.10581, 10.5842, 12.005, 17.6989, 31.5715,
     &      57.746, 85.7318, 265.117, 381.817, 64.2087/)
*         print*,'inside else if 1'
*
      elseif((1.05d-3 .ge. met_cut1) .and. (1.05d-3 .le. met_cut2))then
         mBoost = (/9., 12., 19., 30., 40., 55., 100., 150.,250.,560./)
         RZBoost = (/2.90528, 3.4664, 4.54792, 5.89715, 6.90964,
     &     8.23824, 11.4877, 14.4754, 19.4907, 30.4307/)
         RTBoost = (/8.34864, 10.3943, 13.3912, 22.0983, 33.0246,
     &     67.575, 144.839, 152.163, 337.615, 226.371/)
*         print*,'inside else if 2'
*
      elseif( (2.1d-3 .ge. met_cut1) .and. (2.1d-3 .le. met_cut2) )then
         mBoost = (/9., 12., 15., 30., 40., 55., 100., 150.,250.,575./)
         RZBoost = (/3.07199, 3.65601, 4.17698, 6.19243, 7.26206,
     &      8.68025, 12.1746, 15.4243, 21.3749, 37.7439/)
         RTBoost = (/9.12709, 11.3397, 11.9156, 25.9153, 42.9715,
     &      153.085, 109.867, 241.077, 147.582, 195.808/)
*         print*,'inside else if 3'
      else
         print*,'mass and Rtms data for the given metallicity is not
     &      coded for the Boost models.'
      endif
*      print*,'mBoost is ', mBoost
*      print*,'RZBoost is ', RZBoost
*      print*,'RTBoost is ', RTBoost

*
*     We will interpolate in the log-log space.
      rtms = 10**interp(LOG10(m), size(mBoost), LOG10(mBoost),
     &                       LOG10(RTBoost), 'const')
      rzams = 10**interp(LOG10(m), size(mBoost), LOG10(mBoost),
     &                       LOG10(RZBoost), 'const')
      if (rzams .gt. rtms) then
          rtmsBoost = rzams
      else
          rtmsBoost = rtms
      endif
*      print*,'rtmsBoost, rzams, rtms are ', rtmsBoost, rzams, rtms, met
      return
      end
      real*8 function rtmsBpass(m, met)
      implicit none
      real*8 m, met, r0, alpha, met_cut1, met_cut2

      met_cut1 = 0.9*met
      met_cut2 = 1.1*met

      if ( (1.0d-4 .ge. met_cut1) .and. (1.0d-4 .le. met_cut2) ) then
          r0 = 1.43077288d0
          alpha = 6.2151019d-1

      elseif ((2.0d-4 .ge. met_cut1) .and. (2.0d-4 .le. met_cut2)) then
          r0 = 1.43077288d0
          alpha = 6.2151019d-1

      elseif ((1.0d-3 .ge. met_cut1) .and. (1.0d-3 .le. met_cut2)) then
          r0 = 8.1088547d-1
          alpha = 8.2732973d-1

      elseif ((2.0d-3 .ge. met_cut1) .and. (2.0d-3 .le. met_cut2)) then
          r0 = 4.8850513d-1
          alpha = 9.7369477d-1

      else
         print*,'mass and Rtms data for the given metallicity is not
     &      coded for the BPASS models.'
         r0 = 0.0
         alpha = 0.0

      endif

      rtmsBpass = r0*m**alpha

      return
      end function rtmsBpass

*** NEW CODE BETWEEN THESE LINES (END)
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
      real*8 FUNCTION MLalpha(m,lum,kw)
      implicit none
      integer kw
      real*8 m,lum,loggammae,gammae,Xh
*
* A function for calculating the alpha factor for metallicity-dependent
* winds from Grafener et al. 2011
* (MZ 04/08/19)
*

* First, determine the hydrogen fraction as in Graefener et al. 2011
      if(kw.ge.7 .and. kw.lt.10)then
*        WR stars
         Xh = 0.d0
      else
*        Assume metallicities for different kstars
         if(kw.eq.0 .or. kw.eq.1)then
            Xh = 0.7d0
         elseif(kw.eq.2)then
            Xh = 0.6d0
         elseif(kw.eq.3)then
            Xh = 0.5d0
         elseif(kw.eq.4)then
            Xh = 0.4d0
         elseif(kw.eq.5 .or. kw.eq.6)then
            Xh = 0.2d0
        endif
      endif

      loggammae = -4.813d0 + LOG10(1.d0+xh) + LOG10(lum) - LOG10(m)
      gammae = 10.0d0**(loggammae)
      if(gammae.gt.1d0)then
         gammae = 1.d0
      endif

* calculate metallicity scaling parameter alpha
      if(gammae.lt.(2.d0/3.d0))then
         MLalpha = 0.85d0
      elseif(gammae.ge.(2.d0/3.d0))then
         MLalpha = 2.45d0 - (2.4d0 * gammae)
      endif
*
      return
      end
