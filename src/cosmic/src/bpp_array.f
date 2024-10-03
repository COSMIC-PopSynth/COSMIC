***
        SUBROUTINE WRITETAB(jp,tphys,evolve_type,
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
     &                      bin_state,merger_type,tabname)
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

*
* Write results to bpp or bcm array.
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
        CHARACTER*3 tabname
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
            rsunau = 1/aursun
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

* check which table we are writing to and write the appropriate columns
        if (tabname .eq. 'bpp') then
            jp = MIN(900,jp + 1)        ! Why is the 900 limit here??
            do 117, col_ind = 1, n_col_bpp
                bpp(jp,col_ind) = all_cols(col_inds_bpp(col_ind))
117         continue
        else if (tabname .eq. 'bcm') then
            jp = jp + 1
            do 118, col_ind = 1, n_col_bcm
                bcm(jp,col_ind) = all_cols(col_inds_bcm(col_ind))
118         continue
        end if
        END