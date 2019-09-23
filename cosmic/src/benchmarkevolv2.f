        PROGRAM benchmarkevolv2
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

        INTEGER kstar(2)
        REAL*8 z,ecc,tb,tphysf
        REAL*8 mass(2)
        REAL*8 bppout(1000,15)
        REAL*8 bcmout(50000,42)
        REAL*8 dtptmp
        REAL*8 bhspin(2)
        REAL*8 mass0(2),massc(2),menv(2)
        REAL*8 rad(2),epoch(2)
        REAL*8 lumin(2),renv(2),radc(2)
        REAL*8 zpars(20),bkick(20)
        REAL*8 tacc(2),bacc(2),tms(2),B_0(2),ospin(2)
        REAL*8 tphys

        kstar(1) = 1; kstar(2) = 1
        mass(1) = 33.41813720577207
        mass(2) = 27.46995284892487
        tb = 673.3728182337667
        ecc = 0.6402214090190684; z = 0.002; tphysf = 13700
        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.0; lumin = 0.0; massc = 0.0
        radc = 0.0; menv = 0.0; renv = 0.0
        ospin = 0.0; B_0 = 0.0; bacc = 0.0
        tacc = 0.0 ; epoch = 0.0; tms = 0.0
        bhspin = 0.0; tphys = 0.0
        zpars = 0.0; bkick = 0.0

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

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,bppout,bcmout)

        kstar(1) = 1; kstar(2) = 1; mass(1) = 53.4;
        mass(2) = 45.687; tb = 645.353
        ecc = 0.566449715; z = 0.002; tphysf = 9318.775575930065

        mass0(1) = mass(1)
        mass0(2) = mass(2)

        rad = 0.d0; lumin = 0.d0; massc = 0.d0        
        radc = 0.d0; menv = 0.d0; renv = 0.d0
        ospin = 0.d0; B_0 = 0.d0; bacc = 0.d0
        tacc = 0.d0 ; epoch = 0.d0; tms = 0.d0
        bhspin = 0.d0; tphys = 0.d0
        zpars = 0.d0; bkick = 0.d0

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

        CALL evolv2(kstar,mass,tb,ecc,z,tphysf,
     & dtptmp,mass0,rad,lumin,massc,radc,
     & menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     & bhspin,tphys,zpars,bkick,bppout,bcmout)
        END PROGRAM benchmarkevolv2
