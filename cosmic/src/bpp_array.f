***
        SUBROUTINE WRITEBPP(jp,tphys,evolve_type,
     &                      mass1,mass2,kstar1,kstar2,sep,
     &                      tb,ecc,rrl1,rrl2,bkick,
     &                      aj1,aj2,tms1,tms2,
     &                      massc1,massc2,rad1,rad2)
        IMPLICIT NONE
        INCLUDE 'const_bse.h'
*
* Write results to bpp array.
*
*     Author : Scott Coughlin
*     Date :   12th March 2019
*
        REAL*8 bkick(20),mass1,mass2
        REAL*8 evolve_type,sep,tb,ecc,tphys,rrl1,rrl2
        REAL*8 aj1,aj2,tms1,tms2,massc1,massc2,rad1,rad2
        INTEGER jp,jj
        INTEGER kstar1,kstar2
        REAL*8 yeardy
        PARAMETER(yeardy=365.24d0)

        jp = MIN(80,jp + 1)
        bpp(jp,1) = tphys
        bpp(jp,2) = mass1
        bpp(jp,3) = mass2
        bpp(jp,4) = float(kstar1)
        bpp(jp,5) = float(kstar2)
        bpp(jp,6) = sep
        if(tb.le.0.d0)then
* system was disrupted and tb=-1 and should stay that way
            bpp(jp,7) = tb
        else
            bpp(jp,7) = tb*yeardy
        endif
        bpp(jp,8) = ecc
        bpp(jp,9) = rrl1
        bpp(jp,10) = rrl2
        bpp(jp,11) = evolve_type
        bpp(jp,12) = bkick(15)
        bpp(jp,13) = bkick(16)
* when writing the natal kick for one of the objects
* we have four sitations

*       SURVIVES FIRST SN
        if(bkick(1).eq.1.d0.and.bkick(5).eq.0.d0.and.
     &         bkick(9).eq.0.d0)then
*           write the natal kick for primary object
            bpp(jp,14) = bkick(13)
*           record angular change to the orbital angular momentum
            bpp(jp,15) = bkick(18)
*       DISRUPTS FIRST SN
        elseif(bkick(1).eq.1.d0.and.bkick(5).eq.1.d0.and.
     &         bkick(9).eq.0.d0)then
*           write the natal kick for primary object
            bpp(jp,14) = bkick(13)
*       SURVIVES SECOND SN
        elseif(bkick(1).eq.1.d0.and.bkick(5).eq.2.d0.and.
     &         bkick(9).eq.0.d0)then
*           write the natal kick for secondary object
            bpp(jp,14) = bkick(14)
*           record angular change to the orbital angular momentum
            bpp(jp,15) = bkick(19)
*       DISRUPTS SECOND SN
        elseif(bkick(1).eq.1.d0.and.bkick(5).eq.2.d0.and.
     &         bkick(9).eq.2.d0)then
*           write the natal kick for secondary object
            bpp(jp,14) = bkick(14)
*       SECOND SN AFTER SYSTEM DISRUPTION FROM FIRST SN
        elseif(bkick(1).eq.1.d0.and.bkick(5).eq.1.d0.and.
     &         bkick(9).eq.2.d0)then
*           write the natal kick for secondary object
            bpp(jp,14) = bkick(14)

        endif
        DO jj = 13,20
            bkick(jj) = 0.0
        ENDDO
        bpp(jp,16) = aj1
        bpp(jp,17) = aj2
        bpp(jp,18) = tms1
        bpp(jp,19) = tms2
        bpp(jp,20) = massc1
        bpp(jp,21) = massc2
        bpp(jp,22) = rad1
        bpp(jp,23) = rad2
        END

***
        SUBROUTINE WRITEBCM(ip,tphys,kstar_1,mass0_1,mass_1,
     &                      lumin_1,rad_1,teff_1,massc_1,
     &                      radc_1,menv_1,renv_1,epoch_1,
     &                      ospin_1,deltam_1,RROL_1,kstar_2,mass0_2,
     &                      mass_2,lumin_2,rad_2,teff_2,massc_2,radc_2,
     &                      menv_2,renv_2,epoch_2,ospin_2,deltam_2,
     &                      RROL_2,porb,sep,ecc,B_0_1,B_0_2,SNkick_1,
     &                      SNkick_2, Vsys_final,SNtheta_final,
     &                      SN_1,SN_2,bin_state,merger_type)
        IMPLICIT NONE
        INCLUDE 'const_bse.h'
*
* Write results to bcm array.
*
*     Author : Scott Coughlin
*     Date :   12th March 2019
*
        REAL*8 tphys,mass0_1,mass_1,lumin_1,rad_1,teff_1
        REAL*8 massc_1,radc_1,menv_1,renv_1,epoch_1
        REAL*8 ospin_1,deltam_1,RROL_1
        REAL*8 mass0_2,mass_2,lumin_2,rad_2,teff_2,massc_2
        REAL*8 radc_2,menv_2,renv_2,epoch_2,ospin_2,deltam_2
        REAL*8 RROL_2,porb,sep,ecc,B_0_1,B_0_2
        REAL*8 SNkick_1,SNkick_2,Vsys_final,SNtheta_final
        INTEGER kstar_1,kstar_2,SN_1,SN_2,bin_state,merger_type
        INTEGER ip
        REAL*8 yeardy
        PARAMETER(yeardy=365.24d0)

        ip = ip + 1
        bcm(ip,1) = tphys
        bcm(ip,2) = float(kstar_1)
        bcm(ip,3) = mass0_1
        bcm(ip,4) = mass_1
        bcm(ip,5) = lumin_1
        bcm(ip,6) = rad_1
        bcm(ip,7) = teff_1
        bcm(ip,8) = massc_1
        bcm(ip,9) = radc_1
        bcm(ip,10) = menv_1
        bcm(ip,11) = renv_1
        bcm(ip,12) = epoch_1
        bcm(ip,13) = ospin_1
        bcm(ip,14) = deltam_1
        bcm(ip,15) = RROL_1
        bcm(ip,16) = float(kstar_2)
        bcm(ip,17) = mass0_2
        bcm(ip,18) = mass_2
        bcm(ip,19) = lumin_2
        bcm(ip,20) = rad_2
        bcm(ip,21) = teff_2
        bcm(ip,22) = massc_2
        bcm(ip,23) = radc_2
        bcm(ip,24) = menv_2
        bcm(ip,25) = renv_2
        bcm(ip,26) = epoch_2
        bcm(ip,27) = ospin_2
        bcm(ip,28) = deltam_2
        bcm(ip,29) = RROL_2
        if(porb.le.0.d0)then
* system was disrupted and porb=-1 and should stay that way
            bcm(ip,30) = porb
        else
            bcm(ip,30) = porb*yeardy
        endif
        bcm(ip,31) = sep
        bcm(ip,32) = ecc
        bcm(ip,33) = B_0_1
        bcm(ip,34) = B_0_2
        bcm(ip,35) = SNkick_1
        bcm(ip,36) = SNkick_2
        if(bin_state.ne.2.d0)then
            bcm(ip,37) = Vsys_final
            bcm(ip,38) = SNtheta_final
        endif
        bcm(ip,39) = float(SN_1)
        bcm(ip,40) = float(SN_2)
        bcm(ip,41) = bin_state
        bcm(ip,42) = merger_type

        END
