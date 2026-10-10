      SUBROUTINE zcnsts(z,zpars)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 z,zpars(20)
      integer :: ierr

      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_zcnsts',using_METISSE
*         METISSE (re)processes the tracks when zpars are all zero,
*         so first reset the controls it would read from its input files
          if (all(zpars.eq.0.d0)) call reset_metisse_controls()
          CALL METISSE_zcnsts(z,zpars,ierr)
           if (ierr/=0) call assign_error()
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_zcnsts'
          CALL SSE_zcnsts(z,zpars)
      endif

      END
