      SUBROUTINE zcnsts(z,zpars)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 z,zpars(20)
      
      if (using_METISSE) then
          !WRITE(*,*) 'Calling METISSE_zcnsts',using_METISSE
          !SSE_zcnsts also sets some coefficients used in gntage
          !updating gntage for METISSE will remove the need to call this
          CALL SSE_zcnsts(z,zpars)
          CALL METISSE_zcnsts(z,zpars)
          
      elseif (using_SSE) then
          !WRITE(*,*) 'Calling SSE_zcnsts'
          CALL SSE_zcnsts(z,zpars)
      endif

      END
