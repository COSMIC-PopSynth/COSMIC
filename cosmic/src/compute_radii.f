***
      SUBROUTINE compute_r(mass,z,rad,num)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      CALL zcnsts(z,zpars)
      
      do 10 , k = 1,num
         age = tphys - epoch(k)
         mc = 0.d0 
         lumin(1) = 1.0d-10
         rad = 1.0d-10
         massc = 0.d0
         dmt = 0.d0
         dmr = 0.d0
         rc = radc(k)
         CALL star(kstar,mass0,mass(k),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(k),age,mass(k),tm,tn,tscls,lums,GB,zpars,
     &               rad,lum,kstar,mc,rc,me,re,k2,bhspin,k)

  10 continue
