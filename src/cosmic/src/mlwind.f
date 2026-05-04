      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      integer kw,id
      real*8 lum,r,mt,mc,rl,z
      
      real*8 SSE_mlwind, METISSE_mlwind
      external SSE_mlwind, METISSE_mlwind
    
      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_mlwind'
          mlwind = METISSE_mlwind(kw,lum,r,mt,mc,rl,z,id)
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_mlwind'
          mlwind = SSE_mlwind(kw,lum,r,mt,mc,rl,z)
      endif

      END
