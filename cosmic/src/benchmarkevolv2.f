        PROGRAM benchmarkevolv2
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

        INTEGER kstar1,kstar2
        REAL*8 z,ecc,tb,tphysf,mass1,mass2
        REAL*8 bppout(1000,15)
        REAL*8 bcmout(50000,42)
        REAL*8 dtptmp
        REAL*8 bhspin(2)

        kstar1 = 1; kstar2 = 1; mass1 = 33.41813720577207;
        mass2 = 27.46995284892487; tb = 673.3728182337667
        ecc = 0.6402214090190684; z = 0.002; tphysf = 13700
        bhspin = 0.d0
        neta = 0.5; bwind = 0.0; hewind = 1.0
        alpha1 = 1.0; lambdaf = 1.0; ceflag = 0
        tflag = 1; ifflag = 0; wdflag = 0
        pisn = 45.0; bhflag = 0; nsflag = 3
        cekickflag = 0; cemergeflag = 0; cehestarflag = 0
        mxns = 3.0; pts1 = 0.05; pts2 = 0.01; pts3 = 0.02
        ecsn = 2.5; ecsn_mlow = 1.6; aic = 1; ussn = 0
        sigma = 265.0; sigmadiv = -20.0
        bhsigmafrac = 1.0; polar_kick_angle = 90.0
        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        fprimc_array = 2.d0/21.d0;
        beta = -1.0; xi = 0.5; acc2 = 1.5; epsnov = 0.001
        eddfac = 1.0; gamma = -2.0
        bconst = -3000; CK = -1000; windflag = 3; qcflag = 1
        eddlimflag = 0; dtptmp = 13700.d0; idum1 = 113271
        bhspinflag = 0; bhspinmag=0.d0
        bppout = 0.d0; bcmout = 0.d0

        CALL evolv2(kstar1,kstar2,mass1,mass2,tb,ecc,z,tphysf,
     & dtptmp,bhspin,bppout,bcmout)

        kstar1 = 1; kstar2 = 1; mass1 = 53.4;
        mass2 = 45.687; tb = 645.353
        ecc = 0.566449715; z = 0.002; tphysf = 9318.775575930065
        neta = 0.5; bwind = 0.0; hewind = 1.0
        alpha1 = 1.0; lambdaf = 1.0; ceflag = 0
        tflag = 1; ifflag = 0; wdflag = 0
        pisn = 45.0; bhflag = 0; nsflag = 3
        cekickflag = 0; cemergeflag = 0; cehestarflag = 0
        mxns = 3.0; pts1 = 0.001; pts2 = 0.01; pts3 = 0.02
        ecsn = 2.5; ecsn_mlow = 1.4; aic = 1; ussn = 0
        sigma = 265.0; sigmadiv = -20.0
        bhsigmafrac = 1.0; polar_kick_angle = 90.0
        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        beta = -1.0; xi = 0.5; acc2 = 1.5; epsnov = 0.001
        eddfac = 1.0; gamma = -2.0
        bconst = -3000; CK = -1000; windflag = 3; qcflag = 1
        eddlimflag = 0; dtptmp = 13700.d0; idum1 = 121025
        bppout = 0.d0; bcmout = 0.d0

        CALL evolv2(kstar1,kstar2,mass1,mass2,tb,ecc,z,tphysf,
     & dtptmp,bhspin,bppout,bcmout)
        END PROGRAM benchmarkevolv2
