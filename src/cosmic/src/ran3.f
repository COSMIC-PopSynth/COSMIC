***
      REAL FUNCTION ran3(idum)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* Random number generator from Numerical Recipes, Press et al. pg 272.
*
      INTEGER idum
      INTEGER j,k,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      PARAMETER(im1=2147483563,im2=2147483399,ia1=40014,ia2=40692)
      PARAMETER(iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32)
      DATA idum2/123456789/, iy/0/, ir/ntab*0/
      REAL am
*
      am = 1.0/float(im1)
      imm1 = im1 - 1
      ndiv = 1 + imm1/ntab
*
      if(idum.le.0)then
         idum = MAX(-idum,1)
         idum2 = idum
         do 11 , j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum = idum + im1
            if(j.le.ntab) ir(j) = idum
 11      continue
         iy = ir(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = ir(j) - idum2
      ir(j) = idum
      if(iy.lt.1) iy = iy + imm1
      ran3 = am*iy
*
      RETURN
      END
***



      SUBROUTINE RandomNormal(mean, sigma, idum, result)
* Generate a normally distributed random number with given mean and sigma
* using the Box-Muller transform

      real*8 mean, sigma, result
      integer idum
      real*8 u1, u2, Z0
      
      u1 = ran3(idum)
      u2 = ran3(idum)
      Z0 = SQRT(-2.d0*LOG(u1))*COS(2.d0*3.141592653589793d0*u2)
      result = Z0 * sigma + mean

      RETURN
      END


      SUBROUTINE RandomTruncatedNormal(mu, sigma, idum, lower, upper, x)
* Generate a random number from a truncated normal distribution
* with mean mu, standard deviation sigma, truncated to [lower, upper]
      IMPLICIT NONE
      REAL*8 mu, sigma, lower, upper, x
      INTEGER idum
      REAL ran3
      REAL*8 u, pL, pU, q, z, tiny
      integer max_attempts, attempt
      tiny = 1d-16

      attempt = 0
      max_attempts = 1000

      do
          attempt = attempt + 1
          call RandomNormal(mu, sigma, idum, x)  ! for debugging
          if (x .GE. lower .AND. x .LE. upper)then
             exit
          elseif (attempt .GE. max_attempts) then
            ! use the midpoint if we exceed max attempts
            x = 0.5d0 * (lower + upper)
            exit
          endif
      end do
      RETURN
      END

!     !   WRITE(*,*) ' RandomTruncatedNormal called with: '
!     !   WRITE(*,*) ' mu = ', mu
!     !   WRITE(*,*) ' sigma = ', sigma
!     !   WRITE(*,*) ' lower = ', lower
!     !   WRITE(*,*) ' upper = ', upper

! * map bounds through standard normal cdf
!       call norm_cdf((lower - mu) / sigma, pL)
!       call norm_cdf((upper - mu) / sigma, pU)

!     !   WRITE(*,*) pL, pU

! * clamp in case of extreme tails
!       IF (pL .LT. 0d0) pL = 0d0
!       IF (pU .GT. 1d0) pU = 1d0

! * if interval is numerically tiny in cdf-space, return lower (or midpoint)
!       IF (pU - pL .LE. tiny) THEN
!          x = lower
!          RETURN
!       ENDIF

! * sample uniformly on [pL, pU] then invert
!       u = ran3(idum)
!       q = pL + u*(pU - pL)

!     !   WRITE(*,*) ' sampled q = ', q

! * numerical safety for inverse cdf
!       IF (q .LE. tiny) q = tiny
!       IF (q .GE. 1d0 - tiny) q = 1d0 - tiny

!       call norm_inv(q, z)

!       WRITE(*,*) 'q = ', q, ' z = ', z

!     !   test is PPND7 gives the same result

!       pL = ppnd16(z)
!       write(*,*) ' PPND7(z) = ', pL

!     !   WRITE(*,*) ' inverted z = ', z

!       x = mu + sigma * z

!     !   WRITE(*,*) ' final x = ', x

!       RETURN
!       END


!       SUBROUTINE norm_cdf(z, p)
! *--------------------------------------------------------------
! *  standard normal cdf Phi(z) in real*8
! *  input : z  (real*8)
! *  output: p = Phi(z) (real*8)
! *  method: rational approximation with symmetry; good accuracy
! *          across the range, stable in the tails
! *--------------------------------------------------------------
!       IMPLICIT NONE
!       REAL*8 z, p
!       REAL*8 x, ax, t, phi
!       REAL*8 b1, b2, b3, b4, b5, pp
!       REAL*8 c, tiny

! * coefficients for the as 7.1.26-like approximation
!       b1 = 0.319381530d0
!       b2 = -0.356563782d0
!       b3 = 1.781477937d0
!       b4 = -1.821255978d0
!       b5 = 1.330274429d0
!       c  = 0.2316419d0
!       tiny = 1d-300

!       x  = z
!       ax = ABS(x)

! * fast exits to avoid underflow/overflow in extreme tails
!       IF (x .GE. 8d0) THEN
!          p = 1d0
!          RETURN
!       ELSEIF (x .LE. -8d0) THEN
!          p = 0d0
!          RETURN
!       ENDIF

! * gaussian pdf at |x|
!       phi = 0.39894228040143267794d0 * EXP(-0.5d0*ax*ax)

! * core approximation for upper tail at +|x|
!       t = 1d0 / (1d0 + c*ax)
!       pp = 1d0 - phi * (((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t)

! * reflect for negative x
!       IF (x .LT. 0d0) THEN
!          p = 1d0 - pp
!       ELSE
!          p = pp
!       ENDIF

! * clamp to [0,1] defensively
!       IF (p .LT. 0d0) p = 0d0
!       IF (p .GT. 1d0) p = 1d0
!       IF (p .LT. tiny) p = 0d0
!       RETURN
!       END


      SUBROUTINE norm_inv(p, z)
*--------------------------------------------------------------
*  Compute the inverse CDF (quantile) of the standard normal.
*  Input :  p  in (0,1)
*  Output:  z = Phi^{-1}(p)
*  Algorithm: Peter Acklam's rational approximation
*  Accuracy:  ~1d-9 in REAL*8 precision
*--------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 p, z, q, r
      REAL*8 a1,a2,a3,a4,a5,a6
      REAL*8 b1,b2,b3,b4,b5
      REAL*8 c1,c2,c3,c4,c5,c6
      REAL*8 d1,d2,d3,d4
      REAL*8 tiny

      tiny = 1d-16
      IF (p .LE. 0d0) THEN
         z = -1d300
         RETURN
      ELSEIF (p .GE. 1d0) THEN
         z = 1d300
         RETURN
      ENDIF

* coefficients
      a1 = -3.969683028665376d0
      a2 =  2.209460984245205d0
      a3 = -2.759285104469687d0
      a4 =  1.383577518672690d0
      a5 = -3.066479806614716d-1
      a6 =  2.506628277459239d-2

      b1 = -5.447609879822406d0
      b2 =  1.615858368580409d1
      b3 = -1.556989798598866d1
      b4 =  6.680131188771972d0
      b5 = -1.328068155288572d0

      c1 = -7.784894002430293d-3
      c2 = -3.223964580411365d-1
      c3 = -2.400758277161838d0
      c4 = -2.549732539343734d0
      c5 =  4.374664141464968d0
      c6 =  2.938163982698783d0

      d1 =  7.784695709041462d-3
      d2 =  3.224671290700398d-1
      d3 =  2.445134137142996d0
      d4 =  3.754408661907416d0

* main rational approximation
      IF (p .LT. 0.02425d0) THEN
* lower tail
         q = SQRT(-2d0*LOG(p))
         z = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
     &       ((((d1*q + d2)*q + d3)*q + d4)*q + 1d0)
         z = -z

      ELSEIF (p .GT. 1d0 - 0.02425d0) THEN
* upper tail
         q = SQRT(-2d0*LOG(1d0 - p))
         z = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
     &       ((((d1*q + d2)*q + d3)*q + d4)*q + 1d0)

      ELSE
* central region
         q = p - 0.5d0
         r = q*q
         z = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q /
     &       (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1d0)
      ENDIF

      RETURN
      END

      function ppnd16 ( p, ifault )

c*********************************************************************72
c
cc PPND16 produces the normal deviate value corresponding to lower tail area = P.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**16.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    Michael Wichura
c
c  Reference:
c
c    Michael Wichura,
c    Algorithm AS 241:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 37, Number 3, 1988, pages 477-484.
c
c  Parameters:
c
c    Input, double precision P, the value of the cumulative probability 
c    densitity function.  0 < P < 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, P <= 0 or P >= 1.
c
c    Output, double precision PPND16, the normal deviate value with the 
c    property that the probability of a standard normal deviate being 
c    less than or equal to PPND16 is P.
c
      implicit none

      double precision a0
      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision a6
      double precision a7
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision c5
      double precision c6
      double precision c7
      double precision const1
      double precision const2
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      double precision d5
      double precision d6
      double precision d7
      double precision e0
      double precision e1
      double precision e2
      double precision e3
      double precision e4
      double precision e5
      double precision e6
      double precision e7
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision f5
      double precision f6
      double precision f7
      integer ifault
      double precision p
      double precision ppnd16
      double precision q
      double precision r
      double precision split1
      double precision split2

      parameter ( a0 = 3.3871328727963666080D+00 )
      parameter ( a1 = 1.3314166789178437745D+02 )
      parameter ( a2 = 1.9715909503065514427D+03 )
      parameter ( a3 = 1.3731693765509461125D+04 )
      parameter ( a4 = 4.5921953931549871457D+04 )
      parameter ( a5 = 6.7265770927008700853D+04 )
      parameter ( a6 = 3.3430575583588128105D+04 )
      parameter ( a7 = 2.5090809287301226727D+03 )
      parameter ( b1 = 4.2313330701600911252D+01 )
      parameter ( b2 = 6.8718700749205790830D+02 )
      parameter ( b3 = 5.3941960214247511077D+03 )
      parameter ( b4 = 2.1213794301586595867D+04 )
      parameter ( b5 = 3.9307895800092710610D+04 )
      parameter ( b6 = 2.8729085735721942674D+04 )
      parameter ( b7 = 5.2264952788528545610D+03 )
      parameter ( c0 = 1.42343711074968357734D+00 )
      parameter ( c1 = 4.63033784615654529590D+00 )
      parameter ( c2 = 5.76949722146069140550D+00 )
      parameter ( c3 = 3.64784832476320460504D+00 )
      parameter ( c4 = 1.27045825245236838258D+00 )
      parameter ( c5 = 2.41780725177450611770D-01 )
      parameter ( c6 = 2.27238449892691845833D-02 )
      parameter ( c7 = 7.74545014278341407640D-04 )
      parameter ( const1 = 0.180625D+00 )
      parameter ( const2 = 1.6D+00 )
      parameter ( d1 = 2.05319162663775882187D+00 )
      parameter ( d2 = 1.67638483018380384940D+00 )
      parameter ( d3 = 6.89767334985100004550D-01 )
      parameter ( d4 = 1.48103976427480074590D-01 )
      parameter ( d5 = 1.51986665636164571966D-02 )
      parameter ( d6 = 5.47593808499534494600D-04 )
      parameter ( d7 = 1.05075007164441684324D-09 )
      parameter ( e0 = 6.65790464350110377720D+00 )
      parameter ( e1 = 5.46378491116411436990D+00 )
      parameter ( e2 = 1.78482653991729133580D+00 )
      parameter ( e3 = 2.96560571828504891230D-01 )
      parameter ( e4 = 2.65321895265761230930D-02 )
      parameter ( e5 = 1.24266094738807843860D-03 )
      parameter ( e6 = 2.71155556874348757815D-05 )
      parameter ( e7 = 2.01033439929228813265D-07 )
      parameter ( f1 = 5.99832206555887937690D-01 )
      parameter ( f2 = 1.36929880922735805310D-01 )
      parameter ( f3 = 1.48753612908506148525D-02 )
      parameter ( f4 = 7.86869131145613259100D-04 )
      parameter ( f5 = 1.84631831751005468180D-05 )
      parameter ( f6 = 1.42151175831644588870D-07 )
      parameter ( f7 = 2.04426310338993978564D-15 )
      parameter ( split1 = 0.425D+00 )
      parameter ( split2 = 5.D+00 )

      ifault = 0
      q = p - 0.5D+00

      if ( dabs ( q ) .le. split1 ) then

        r = const1 - q * q

        ppnd16 = q * (((((((
     &      a7   * r 
     &    + a6 ) * r 
     &    + a5 ) * r 
     &    + a4 ) * r 
     &    + a3 ) * r 
     &    + a2 ) * r 
     &    + a1 ) * r 
     &    + a0 ) / (((((((
     &      b7   * r 
     &    + b6 ) * r 
     &    + b5 ) * r 
     &    + b4 ) * r 
     &    + b3 ) * r 
     &    + b2 ) * r 
     &    + b1 ) * r 
     &    + 1.0D+00 )

      else

        if ( q .lt. 0.0D+00 ) then
          r = p
        else
          r = 1.0D+00 - p
        end if

        if ( r .le. 0.0D+00 ) then
          ifault = 1
          ppnd16 = 0.0D+00
          return
        end if

        r = dsqrt ( - dlog ( r ) )

        if ( r .le. split2 ) then

          r = r - const2

          ppnd16 = (((((((
     &        c7   * r 
     &      + c6 ) * r 
     &      + c5 ) * r 
     &      + c4 ) * r 
     &      + c3 ) * r 
     &      + c2 ) * r 
     &      + c1 ) * r 
     &      + c0 ) / (((((((
     &        d7   * r 
     &      + d6 ) * r 
     &      + d5 ) * r 
     &      + d4 ) * r 
     &      + d3 ) * r 
     &      + d2 ) * r 
     &      + d1 ) * r 
     &      + 1.0D+00 )

        else

          r = r - split2

          ppnd16 = (((((((
     &        e7   * r 
     &      + e6 ) * r 
     &      + e5 ) * r 
     &      + e4 ) * r 
     &      + e3 ) * r 
     &      + e2 ) * r 
     &      + e1 ) * r 
     &      + e0 ) / (((((((
     &        f7   * r 
     &      + f6 ) * r 
     &      + f5 ) * r 
     &      + f4 ) * r 
     &      + f3 ) * r 
     &      + f2 ) * r 
     &      + f1 ) * r 
     &      + 1.0D+00 )

        end if

        if ( q .lt. 0.0D+00 ) then
          ppnd16 = - ppnd16
        end if

      end if

      return
      end