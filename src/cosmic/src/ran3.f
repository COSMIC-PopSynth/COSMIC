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



      SUBROUTINE RandomNormal(mean, sigma, idum, lower, upper, result)
* Generate a normally distributed random number with given mean and sigma
* using the Box-Muller transform, redrawing if the number is outside
* the specified lower and upper bounds.

      real*8 mean, sigma, result, lower, upper
      integer idum
      real*8 u1, u2, Z0

      do
         u1 = ran3(idum)
         u2 = ran3(idum)
         Z0 = SQRT(-2.d0*LOG(u1))*COS(2.d0*3.141592653589793d0*u2)
         result = Z0 * sigma + mean
         if(result.ge.lower .and. result.le.upper) exit
      enddo

      RETURN
      END
