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


      SUBROUTINE RandomTruncatedNormal(mu, sigma, idum, lower, upper,
     &                                 n_try, x)
* Draw from a normal distribution with mean mu and standard deviation
* sigma truncated to [lower, upper], by rejection sampling.
*
* If n_try > 0, this always makes exactly 2 * n_try calls to RAN3 (2 per
* candidate, with the calls for any candidates left after one is
* accepted still made). This keeps the random number stream the same
* whatever the outcome, so that later draws don't depend on how many
* candidates were needed. In the (vanishingly unlikely) case that all
* n_try candidates fall outside [lower, upper], the last one is
* clamped to the range.
*
* If n_try <= 0, candidates are drawn until one is accepted (so the
* number of calls to RAN3 varies), up to a maximum of 1000, after which
* the midpoint of the range is used.
      IMPLICIT NONE
      REAL*8 mu, sigma, lower, upper, x, y
      INTEGER idum, n_try, attempt, max_attempts
      PARAMETER(max_attempts=1000)
      LOGICAL found
      REAL ran3, pad
      EXTERNAL ran3

      if(n_try.gt.0)then
         found = .false.
         y = mu
         do attempt = 1, n_try
            if(found)then
               pad = ran3(idum)
               pad = ran3(idum)
            else
               call RandomNormal(mu, sigma, idum, y)
               if(y.ge.lower.and.y.le.upper)then
                  x = y
                  found = .true.
               endif
            endif
         end do
         if(.not.found) x = MIN(MAX(y, lower), upper)
      else
         attempt = 0
         do
            attempt = attempt + 1
            call RandomNormal(mu, sigma, idum, x)
            if(x.ge.lower.and.x.le.upper)then
               exit
            elseif(attempt.ge.max_attempts)then
               x = 0.5d0 * (lower + upper)
               exit
            endif
         end do
      endif

      RETURN
      END
