      SUBROUTINE zcnsts(z,zpars,path_to_tracks,path_to_he_tracks)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 z,zpars(20)
      CHARACTER*256 path_to_tracks,path_to_he_tracks
      integer :: ierr

      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_zcnsts',using_METISSE
          CALL METISSE_zcnsts(z,zpars,path_to_tracks,
     &     path_to_he_tracks,ierr)
           if (ierr/=0) call assign_error()
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_zcnsts'
          CALL SSE_zcnsts(z,zpars)
      endif

      END
