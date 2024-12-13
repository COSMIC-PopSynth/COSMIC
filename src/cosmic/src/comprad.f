***
      SUBROUTINE compute_r(mass,z,num,rad,
     &           path_to_tracks,path_to_he_tracks)
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
      real*8 mt,tm,tn,mass0,age,lum,mc,rc,me,re,dtm
      REAL*8 tscls(20),lums(10),GB(10),zpars(20),k2,bhspin
      CHARACTER*256 path_to_tracks,path_to_he_tracks

***
* f2py directives go here; we'll return the radii as a 10^5 array
***

Cf2py intent(in) mass
Cf2py intent(in) z
Cf2py intent(in) num
Cf2py intent(out) rad

      if(using_METISSE.eq.1) CALL initialize_front_end('cosmic')
      CALL zcnsts(z,zpars,path_to_tracks,path_to_he_tracks)
      
      if(using_METISSE.eq.1) call allocate_track(num,mass)

      
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
         dtm = 0.d0
         CALL star(kstar,mass0,mt,tm,tn,tscls,lums,GB,zpars,dtm,k)
         CALL hrdiag(mass0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rad(k),lum,kstar,mc,rc,me,re,k2,bhspin,k)
         
  10  continue
  
        
        if (using_METISSE.eq.1) call dealloc_track()

      END SUBROUTINE compute_r
