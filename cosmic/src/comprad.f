***
      SUBROUTINE compute_r(mass,z,num,rad)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

***
* In Fortran, we can't dynamically allocate arrays
* Instead, we declare a stanard array with 10^5 elements
* and pass everything from cosmic 10^5 at a time
*
* for safety, we'll pass the number of stars, num, and stop the
* for loops there
***


      integer loop
      PARAMETER(loop=100000)
      real*8 mass(loop),rad(loop),z
      integer k,kstar,num
      real*8 mt,tm,tn,mass0,age,lum,mc,rc,me,re
      REAL*8 tscls(20),lums(10),GB(10),zpars(20),k2,bhspin

***
* f2py directives go here; we'll return the radii as a 10^5 array
***

Cf2py intent(in) mass
Cf2py intent(in) z
Cf2py intent(in) num
Cf2py intent(out) rad

      CALL zcnsts(z,zpars)
      
***
* Then just loop through everything 
***
      do 10 , k = 1,num
         age = 0.0
         mc = 0.d0 
         tm = 0.d0
         me = 0.d0 
         re = 0.d0 
         lum = 1.0d-10
         mass0 = mass(k) 
         mt = mass(k)
         kstar = 0
         if(mt.ge.0.7) kstar = 1
         bhspin = 0.d0
         rc = 0.d0 
         CALL star(kstar,mass0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rad(k),lum,kstar,mc,rc,me,re,k2,bhspin,k)

  10  continue

      END SUBROUTINE compute_r
