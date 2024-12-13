      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  bhspin,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      integer kw,id
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 bhspin
      real*8 r,lum,mc,rc,menv,renv,k2,mcx
      
      if (using_METISSE.eq.1) then
          !WRITE(*,*) 'Calling METISSE_hrdiag'
          CALL METISSE_hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  mcx,id)
           ! get_bhspin is defined in assign_commons_cosmic.f90
           if (kw==14) CALL get_bhspin(bhspin,id)
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_hrdiag'
          CALL SSE_hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  bhspin,id)
      endif

      END
