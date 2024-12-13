      SUBROUTINE deltat(kw,age,tm,tn,tscls,dt,dtr,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      INTEGER kw,id
      REAL*8 age,tm,tn,tscls(20)
      REAL*8 dt,dtr
      
      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_deltat'
          CALL METISSE_deltat(id,age,dt,dtr)
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_deltat'
          CALL SSE_deltat(kw,age,tm,tn,tscls,dt,dtr)
      endif

      END
