        PROGRAM benchmarkevolv2
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

        INTEGER kstar(2)
        REAL*8 z,ecc,tb,tphysf
        REAL*8 mass(2)
        REAL*8 dtptmp
        REAL*8 bhspin(2)
        REAL*8 mass0(2),massc(2),menv(2)
        REAL*8 rad(2),epoch(2)
        REAL*8 lumin(2),renv(2),radc(2)
        REAL*8 zpars(20),kick_info(2,17)
        REAL*8 tacc(2),bacc(2),tms(2),B_0(2),ospin(2),bkick(20)
        REAL*8 tphys
        REAL*8 kick_info_out(2,17)
        INTEGER bpp_index_out,bcm_index_out


        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 0.0; kstar(2) = 0.0
        mass(1) = 0.5
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 11.555555555555555
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 11.555555555555555
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 11.555555555555555
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 11.555555555555555
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 11.555555555555555
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 3.263888888888889
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 3.263888888888889
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 3.263888888888889
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 3.263888888888889
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 6.027777777777778
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 6.027777777777778
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 6.027777777777778
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 6.027777777777778
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 6.027777777777778
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 8.791666666666666
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 8.791666666666666
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 8.791666666666666
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 8.791666666666666
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 8.791666666666666
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 11.555555555555555
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 11.555555555555555
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 11.555555555555555
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 11.555555555555555
        mass(2) = 11.555555555555555
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 22.61111111111111
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 22.61111111111111
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 22.61111111111111
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 22.61111111111111
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 22.61111111111111
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 6.027777777777778
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 6.027777777777778
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 6.027777777777778
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 6.027777777777778
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 6.027777777777778
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 11.555555555555555
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 11.555555555555555
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 11.555555555555555
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 11.555555555555555
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 11.555555555555555
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 17.083333333333332
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 17.083333333333332
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 17.083333333333332
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 17.083333333333332
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 17.083333333333332
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 22.61111111111111
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 22.61111111111111
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 22.61111111111111
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 22.61111111111111
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 22.61111111111111
        mass(2) = 22.61111111111111
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 33.666666666666664
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 33.666666666666664
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 33.666666666666664
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 33.666666666666664
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 33.666666666666664
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 8.791666666666666
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 8.791666666666666
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 8.791666666666666
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 8.791666666666666
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 8.791666666666666
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 17.083333333333332
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 17.083333333333332
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 17.083333333333332
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 17.083333333333332
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 17.083333333333332
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 25.375
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 25.375
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 25.375
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 25.375
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 25.375
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 33.666666666666664
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 33.666666666666664
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 33.666666666666664
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 33.666666666666664
        mass(2) = 33.666666666666664
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 44.72222222222222
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 44.72222222222222
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 44.72222222222222
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 44.72222222222222
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 44.72222222222222
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 11.555555555555555
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 11.555555555555555
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 11.555555555555555
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 11.555555555555555
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 11.555555555555555
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 22.61111111111111
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 22.61111111111111
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 22.61111111111111
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 22.61111111111111
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 22.61111111111111
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 33.666666666666664
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 33.666666666666664
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 33.666666666666664
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 33.666666666666664
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 33.666666666666664
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 44.72222222222222
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 44.72222222222222
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 44.72222222222222
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 44.72222222222222
        mass(2) = 44.72222222222222
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 55.77777777777778
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 55.77777777777778
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 55.77777777777778
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 55.77777777777778
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 55.77777777777778
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 14.319444444444445
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 14.319444444444445
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 14.319444444444445
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 14.319444444444445
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 14.319444444444445
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 28.13888888888889
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 28.13888888888889
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 28.13888888888889
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 28.13888888888889
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 28.13888888888889
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 41.958333333333336
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 41.958333333333336
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 41.958333333333336
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 41.958333333333336
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 41.958333333333336
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 55.77777777777778
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 55.77777777777778
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 55.77777777777778
        mass(2) = 55.77777777777778
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 66.83333333333333
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 66.83333333333333
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 66.83333333333333
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 66.83333333333333
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 66.83333333333333
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 17.083333333333332
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 17.083333333333332
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 17.083333333333332
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 17.083333333333332
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 17.083333333333332
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 33.666666666666664
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 33.666666666666664
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 33.666666666666664
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 33.666666666666664
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 33.666666666666664
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 50.25
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 50.25
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 50.25
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 50.25
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 50.25
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 66.83333333333333
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 66.83333333333333
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 66.83333333333333
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 66.83333333333333
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 66.83333333333333
        mass(2) = 66.83333333333333
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 77.88888888888889
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 77.88888888888889
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 77.88888888888889
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 77.88888888888889
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 77.88888888888889
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 19.84722222222222
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 19.84722222222222
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 19.84722222222222
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 19.84722222222222
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 19.84722222222222
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 39.19444444444444
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 39.19444444444444
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 39.19444444444444
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 39.19444444444444
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 39.19444444444444
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 58.541666666666664
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 58.541666666666664
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 58.541666666666664
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 58.541666666666664
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 58.541666666666664
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 77.88888888888889
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 77.88888888888889
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 77.88888888888889
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 77.88888888888889
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 77.88888888888889
        mass(2) = 77.88888888888889
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 88.94444444444444
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 88.94444444444444
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 88.94444444444444
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 88.94444444444444
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 88.94444444444444
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 22.61111111111111
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 22.61111111111111
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 22.61111111111111
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 22.61111111111111
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 22.61111111111111
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 44.72222222222222
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 44.72222222222222
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 44.72222222222222
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 44.72222222222222
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 44.72222222222222
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 66.83333333333333
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 66.83333333333333
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 66.83333333333333
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 66.83333333333333
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 66.83333333333333
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 88.94444444444444
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 88.94444444444444
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 88.94444444444444
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 88.94444444444444
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 88.94444444444444
        mass(2) = 88.94444444444444
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 100.0
        mass(2) = 0.5
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 100.0
        mass(2) = 0.5
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 100.0
        mass(2) = 0.5
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 100.0
        mass(2) = 0.5
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 0.0
        mass(1) = 100.0
        mass(2) = 0.5
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 25.375
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 25.375
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 25.375
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 25.375
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 25.375
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 50.25
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 50.25
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 50.25
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 50.25
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 50.25
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 75.125
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 75.125
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 75.125
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 75.125
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 75.125
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 100.0
        tb = 1000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 100.0
        tb = 3162.2776601683795
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 100.0
        tb = 10000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 100.0
        tb = 31622.776601683792
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

    
        kstar(1) = 1.0; kstar(2) = 1.0
        mass(1) = 100.0
        mass(2) = 100.0
        tb = 100000.0
        ecc = 0.0; z = 0.014; tphysf = 13700.0
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; kick_info = 0.0; bkick = 0.0

        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        natal_kick_array(1, 5) = 0.d0
        natal_kick_array(2, 5) = 0.d0
        fprimc_array = 2.d0/21.d0;

        xi = 1.0
        bhflag = 1
        neta = 0.5
        windflag = 3
        wdflag = 1
        alpha1 = 1.0
        pts1 = 0.001
        pts3 = 0.02
        pts2 = 0.01
        epsnov = 0.001
        hewind = 0.5
        ck = 1000
        bwind = 0.0
        lambdaf = 0.0
        mxns = 2.5
        beta = -1.0
        tflag = 1
        acc2 = 1.5
        remnantflag = 0
        ceflag = 0
        eddfac = 1.0
        ifflag = 0
        bconst = 3000
        sigma = 265.0
        gamma = -2.0
        pisn = 45.0
        bhsigmafrac = 1.0
        polar_kick_angle= 90
        cekickflag = 2
        cehestarflag = 0
        cemergeflag = 0
        ecsn = 2.5
        ecsn_mlow = 1.4
        aic = 1
        ussn = 0
        sigmadiv = -20.0
        qcflag = 4
        eddlimflag = 0
        bhspinflag = 0
        bhspinmag = 0.0
        rejuv_fac = 1.0
        rejuvflag = 0
        htpmb = 1
        ST_cr = 1
        ST_tide = 0
        bdecayfac = 1
        grflag = 1
        rembar_massloss = 0.5
        kickflag = 0
        zsun = 0.014
        bhms_coll_flag = 0
        don_lim = -1
        acc_lim = -1

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,kick_info,
     & bpp_index_out,bcm_index_out,kick_info_out)

            END PROGRAM benchmarkevolv2
