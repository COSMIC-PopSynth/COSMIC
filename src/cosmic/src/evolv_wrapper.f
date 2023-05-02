      
      SUBROUTINE evolv2_global(z,zpars, alpha2, acc_lim1, q3, q4,
     \ natal_kick)
         
         implicit none
         INCLUDE 'const_bse.h'
         integer i,kw,kw2,kstar(2),j,k,time,idum
         INTEGER bcm_index_out, bpp_index_out
   
         REAL*8 mass0(2),mass(2),z,zpars(20),kick_info(2,17)
         REAL*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
         REAL*8 rad(2),lum(2),ospin(2),kick_info_out(2,17)
         REAL*8 massc(2),radc(2),menv(2),renv(2),bhspin(2)
         REAL*8 tb,ecc,yearsc,bkick(20)
         REAL*8 B_0(2),bacc(2),tacc(2),xip,xihold
         REAL*8 alpha2,acc_lim1, q3, q4
         REAL*8 natal_kick(2,5)
   
         PARAMETER(yearsc=3.1557d+07)
         
         dtp = 0.0
         pts1 = 0.01
         pts2 = 0.01
         pts3 = 0.02
         zsun = 0.014
         windflag = 3
         eddlimflag = 0
         neta = 0.5
         bwind = 0.0
         hewind = 0.5
         beta = -1
         xi = 0.5
         acc2 = 1.5
         alpha1 = alpha2
         lambdaf = 0.0
         ceflag = 1
         cekickflag = 2
         cemergeflag = 1
         cehestarflag = 0
         qcflag = 5
         kickflag = 0
         sigma = 0.0
         bhflag = 1
         bhsigmafrac = 1.0
         sigmadiv = -20.0
         ecsn = 2.25
         ecsn_mlow = 1.6
         aic = 1
         ussn = 1
         pisn = -2
         polar_kick_angle = 90.0
         remnantflag = 4
         mxns = 3.0
         rembar_massloss = 0.5
         bhspinflag = 0
         bhspinmag = 0.0
         grflag = 1
         eddfac = 1
         gamma = -2.0
         don_lim = -1
         acc_lim = acc_lim1
         tflag = 1
         ST_tide = 1
         ifflag = 0
         wdflag = 1
         epsnov = 0.001
         bdecayfac = 1
         bconst = 3000
         ck = 1000
         rejuv_fac = 1.0
         rejuvflag = 0
         bhms_coll_flag = 0
         htpmb = 1
         ST_cr = 1
   
         idum = 50
         do j=1,2
            do i=1,12
               kick_info(j,i) = 0.d0
            enddo
            do i=1,5
               natal_kick_array(j,i) = natal_kick(j,i)
            enddo
         enddo
   
         do j=1,16
            qcrit_array(j) = 0.d0
            fprimc_array(j) = 2.0/21.0
         enddo
   
         qcrit_array(3) = q3
         qcrit_array(4) = q4
   
         CALL zcnsts(z,zpars)
         CALL instar
         CALL zcnsts(z,zpars)

         WRITE(*,*)qcrit_array(3),qcrit_array(4)
         WRITE(*,*)natal_kick_array
   
         return
         END
   
