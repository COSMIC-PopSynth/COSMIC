        PROGRAM benchmarkevolv2
        IMPLICIT NONE
        INCLUDE 'const_bse.h'

        INTEGER kstar1,kstar2
        REAL*8 z,ecc,tb,tphysf,mass1,mass2
        REAL*8 bppout(1000,15)
        REAL*8 bcmout(50000,42)

        REAL*8 netatmp,bwindtmp,hewindtmp,alpha1tmp,lambdatmp
        REAL*8 mxnstmp,pts1tmp,pts2tmp,pts3tmp,dtptmp
        REAL*8 sigmatmp,bhsigmafractmp,polar_kick_angletmp,betatmp,xitmp
        REAL*8 ecsntmp,ecsn_mlowtmp,sigmadivtmp
        REAL*8 acc2tmp,epsnovtmp,eddfactmp,gammatmp
        REAL*8 bconsttmp,CKtmp,qc_fixed,qcrit_array(16)
        REAL*8 vk1_bcm,vk2_bcm,vsys_bcm,theta_bcm,natal_kick_array(6)
        INTEGER cekickflagtmp,cemergeflagtmp,cehestarflagtmp,ussntmp
        INTEGER ceflagtmp,tflagtmp,ifflagtmp,nsflagtmp,aictmp,qcflagtmp
        INTEGER wdflagtmp,pisntmp,bhflagtmp,windflagtmp,idumtmp

        kstar1 = 1; kstar2 = 1; mass1 = 33.41813720577207;
        mass2 = 27.46995284892487; tb = 673.3728182337667
        ecc = 0.6402214090190684; z = 0.002; tphysf = 13700
        netatmp = 0.5; bwindtmp = 0.0; hewindtmp = 1.0
        alpha1tmp = 1.0; lambdatmp = 1.0; ceflagtmp = 0
        tflagtmp = 1; ifflagtmp = 0; wdflagtmp = 0
        pisntmp = 45.0; bhflagtmp = 0; nsflagtmp = 3
        cekickflagtmp = 0; cemergeflagtmp = 0; cehestarflagtmp = 0
        mxnstmp = 3.0; pts1tmp = 0.05; pts2tmp = 0.01; pts3tmp = 0.02
        ecsntmp = 2.5; ecsn_mlowtmp = 1.6; aictmp = 1; ussntmp = 0
        sigmatmp = 265.0; sigmadivtmp = -20.0
        bhsigmafractmp = 1.0; polar_kick_angletmp = 90.0
        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        betatmp = -1.0; xitmp = 0.5; acc2tmp = 1.5; epsnovtmp = 0.001
        eddfactmp = 1.0; gammatmp = -2.0
        bconsttmp = -3000; CKtmp = -1000; windflagtmp = 3; qcflagtmp = 1
        dtptmp = 13700.d0; idumtmp = 113271
        bppout = 0.d0; bcmout = 0.d0

        CALL evolv2(kstar1,kstar2,mass1,mass2,tb,ecc,z,tphysf,
     & netatmp,bwindtmp,hewindtmp,alpha1tmp,lambdatmp,
     & ceflagtmp,tflagtmp,ifflagtmp,wdflagtmp,pisntmp,
     & bhflagtmp,nsflagtmp,
     & cekickflagtmp,cemergeflagtmp,cehestarflagtmp,
     & mxnstmp,pts1tmp,pts2tmp,pts3tmp,ecsntmp,ecsn_mlowtmp,aictmp,
     & ussntmp,sigmatmp,sigmadivtmp,bhsigmafractmp,polar_kick_angletmp,
     & natal_kick_array,qcrit_array,betatmp,xitmp,
     & acc2tmp,epsnovtmp,eddfactmp,gammatmp,
     & bconsttmp,CKtmp,windflagtmp,qcflagtmp,
     & dtptmp,idumtmp,bppout,bcmout)

        kstar1 = 1; kstar2 = 1; mass1 = 53.4;
        mass2 = 45.687; tb = 645.353
        ecc = 0.566449715; z = 0.002; tphysf = 9318.775575930065
        netatmp = 0.5; bwindtmp = 0.0; hewindtmp = 1.0
        alpha1tmp = 1.0; lambdatmp = 1.0; ceflagtmp = 0
        tflagtmp = 1; ifflagtmp = 0; wdflagtmp = 0
        pisntmp = 45.0; bhflagtmp = 0; nsflagtmp = 3
        cekickflagtmp = 0; cemergeflagtmp = 0; cehestarflagtmp = 0
        mxnstmp = 3.0; pts1tmp = 0.001; pts2tmp = 0.01; pts3tmp = 0.02
        ecsntmp = 2.5; ecsn_mlowtmp = 1.4; aictmp = 1; ussntmp = 0
        sigmatmp = 265.0; sigmadivtmp = -20.0
        bhsigmafractmp = 1.0; polar_kick_angletmp = 90.0
        natal_kick_array = -100.d0; qcrit_array = 0.d0;
        betatmp = -1.0; xitmp = 0.5; acc2tmp = 1.5; epsnovtmp = 0.001
        eddfactmp = 1.0; gammatmp = -2.0
        bconsttmp = -3000; CKtmp = -1000; windflagtmp = 3; qcflagtmp = 1
        dtptmp = 13700.d0; idumtmp = 121025
        bppout = 0.d0; bcmout = 0.d0

        CALL evolv2(kstar1,kstar2,mass1,mass2,tb,ecc,z,tphysf,
     & netatmp,bwindtmp,hewindtmp,alpha1tmp,lambdatmp,
     & ceflagtmp,tflagtmp,ifflagtmp,wdflagtmp,pisntmp,
     & bhflagtmp,nsflagtmp,
     & cekickflagtmp,cemergeflagtmp,cehestarflagtmp,
     & mxnstmp,pts1tmp,pts2tmp,pts3tmp,ecsntmp,ecsn_mlowtmp,aictmp,
     & ussntmp,sigmatmp,sigmadivtmp,bhsigmafractmp,polar_kick_angletmp,
     & natal_kick_array,qcrit_array,betatmp,xitmp,
     & acc2tmp,epsnovtmp,eddfactmp,gammatmp,
     & bconsttmp,CKtmp,windflagtmp,qcflagtmp,
     & dtptmp,idumtmp,bppout,bcmout)
        END PROGRAM benchmarkevolv2
