***
        SUBROUTINE WRITEBPP(jp,tphys,evolve_type,
     &                      mass1,mass2,kstar1,kstar2,sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj1,aj2,tms1,tms2,
     &                      massc1,massc2,rad1,rad2,
     &                      mass0_1,mass0_2,lumin1,lumin2,
     &                      teff1,teff2,radc1,radc2,menv1,
     &                      menv2,renv1,renv2,ospin1,ospin2,
     &                      b_0_1,b_0_2,bacc1,bacc2,tacc1,tacc2,
     &                      epoch1,epoch2,bhspin1,bhspin2,
     &                      deltam_1,deltam_2,SN_1,SN_2,
     &                      bin_state,merger_type)
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

        ! deltam1_bcm,deltam2_bcm,formation(1),formation(2),binstate,mergertype
*
* Write results to bpp array.
*
*     Author : Scott Coughlin, Tom Wagg
*     Date :   12th March 2019, September 2024
*
        REAL*8 mass1,mass2
        REAL*8 evolve_type,sep,tb,ecc,tphys,rrl1,rrl2
        REAL*8 aj1,aj2,tms1,tms2,massc1,massc2,rad1,rad2
        REAL*8 mass0_1,mass0_2,lumin1,lumin2,radc1,radc2
        REAL*8 menv1,menv2,renv1,renv2,ospin1,ospin2
        REAL*8 b_0_1,b_0_2,bacc1,bacc2,tacc1,tacc2,epoch1,epoch2
        REAL*8 bhspin1,bhspin2,teff1,teff2
        REAL*8 deltam_1,deltam_2
        INTEGER SN_1,SN_2,bin_state,merger_type
        REAL*8 tb_write,sep_cubed
        INTEGER jp, col_ind
        INTEGER kstar1,kstar2
        REAL*8 yeardy,aursun,rsunau
        REAL*8 all_cols(49)
        PARAMETER(yeardy=365.24d0,aursun=214.95d0)

        all_cols(1) = tphys
        all_cols(2) = mass1
        all_cols(3) = mass2
        all_cols(4) = float(kstar1)
        all_cols(5) = float(kstar2)
        all_cols(6) = sep
        if(tb.le.0.d0)then
* system was disrupted and tb=-1 and should stay that way
            all_cols(7) = tb
        else
            sep_cubed = (sep*rsunau)*(sep*rsunau)*(sep*rsunau)
            tb_write = sqrt(sep_cubed/(mass1+mass2))
            all_cols(7) = tb_write*yeardy
        endif
        all_cols(8) = ecc
        all_cols(9) = rrl1
        all_cols(10) = rrl2
        all_cols(11) = evolve_type
        all_cols(12) = aj1
        all_cols(13) = aj2
        all_cols(14) = tms1
        all_cols(15) = tms2
        all_cols(16) = massc1
        all_cols(17) = massc2
        all_cols(18) = rad1
        all_cols(19) = rad2
        all_cols(20) = mass0_1
        all_cols(21) = mass0_2
        all_cols(22) = lumin1
        all_cols(23) = lumin2
        all_cols(24) = teff1
        all_cols(25) = teff2
        all_cols(26) = radc1
        all_cols(27) = radc2
        all_cols(28) = menv1
        all_cols(29) = menv2
        all_cols(30) = renv1
        all_cols(31) = renv2
        all_cols(32) = ospin1
        all_cols(33) = ospin2
        all_cols(34) = b_0_1
        all_cols(35) = b_0_2
        all_cols(36) = bacc1
        all_cols(37) = bacc2
        all_cols(38) = tacc1
        all_cols(39) = tacc2
        all_cols(40) = epoch1
        all_cols(41) = epoch2
        all_cols(42) = bhspin1
        all_cols(43) = bhspin2
        all_cols(44) = deltam_1
        all_cols(45) = deltam_2
        all_cols(46) = float(SN_1)
        all_cols(47) = float(SN_2)
        all_cols(48) = bin_state
        all_cols(49) = merger_type


        rsunau = 1/aursun
        jp = MIN(900,jp + 1)

        do 117, col_ind = 1, n_col_bpp
            bpp(jp,col_ind) = all_cols(col_inds_bpp(col_ind))
117     continue
        END

***
        SUBROUTINE WRITEBCM(ip,tphys,kstar_1,mass0_1,mass_1,
     &                      lumin_1,rad_1,teff_1,massc_1,
     &                      radc_1,menv_1,renv_1,epoch_1,
     &                      ospin_1,deltam_1,RRLO_1,kstar_2,mass0_2,
     &                      mass_2,lumin_2,rad_2,teff_2,massc_2,radc_2,
     &                      menv_2,renv_2,epoch_2,ospin_2,deltam_2,
     &                      RRLO_2,porb,sep,ecc,B_0_1,B_0_2,
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
        REAL*8 ospin_1,deltam_1,RRLO_1,porb_write,sep_cubed
        REAL*8 mass0_2,mass_2,lumin_2,rad_2,teff_2,massc_2
        REAL*8 radc_2,menv_2,renv_2,epoch_2,ospin_2,deltam_2
        REAL*8 RRLO_2,porb,sep,ecc,B_0_1,B_0_2
        INTEGER kstar_1,kstar_2,SN_1,SN_2,bin_state,merger_type
        INTEGER ip
        REAL*8 yeardy,aursun,rsunau
        PARAMETER(yeardy=365.24d0,aursun=214.95d0)

        rsunau = 1/aursun

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
        bcm(ip,15) = RRLO_1
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
        bcm(ip,29) = RRLO_2
        if(porb.le.0.d0)then
* system was disrupted and porb=-1 and should stay that way
            bcm(ip,30) = porb
        else
            sep_cubed = (sep*rsunau)*(sep*rsunau)*(sep*rsunau)
            porb_write = sqrt(sep_cubed/(mass_1+mass_2))
            bcm(ip,30) = porb_write*yeardy
        endif
        bcm(ip,31) = sep
        bcm(ip,32) = ecc
        bcm(ip,33) = B_0_1
        bcm(ip,34) = B_0_2
        bcm(ip,35) = float(SN_1)
        bcm(ip,36) = float(SN_2)
        bcm(ip,37) = bin_state
        bcm(ip,38) = merger_type

        END
