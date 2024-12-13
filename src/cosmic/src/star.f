      SUBROUTINE star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20),dtm
      integer kw ,id
      
      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_star'
          CALL METISSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_star'
          CALL SSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
      endif

      END
