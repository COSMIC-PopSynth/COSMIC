*
* checkstate.h
*
      REAL*8 mass_1_state(15,2),mass_2_state(15,2)
      REAL*8 kstar_1_state(15,2),kstar_2_state(15,2)
      REAL*8 sep_state(15,2),porb_state(15,2),ecc_state(15,2)
      REAL*8 RROL_1_state(15,2),RROL_2_state(15,2),evol_type_state(15,2)
      REAL*8 aj_1_state(15,2),aj_2_state(15,2)
      REAL*8 tms_1_state(15,2),tms_2_state(15,2)
      REAL*8 massc_1_state(15,2),massc_2_state(15,2)
      REAL*8 rad_1_state(15,2),rad_2_state(15,2)
      REAL*8 mass0_1_state(15,2),mass0_2_state(15,2)
      REAL*8 lum_1_state(15,2),lum_2_state(15,2)
      REAL*8 radc_1_state(15,2),radc_2_state(15,2)
      REAL*8 menv_1_state(15,2),menv_2_state(15,2)
      REAL*8 renv_1_state(15,2),renv_2_state(15,2)
      REAL*8 omega_spin_1_state(15,2),omega_spin_2_state(15,2)
      REAL*8 B0_1_state(15,2),B0_2_state(15,2)
      REAL*8 bacc_1_state(15,2),bacc_2_state(15,2)
      REAL*8 tacc_1_state(15,2),tacc_2_state(15,2)
      REAL*8 epoch_1_state(15,2),epoch_2_state(15,2)
      REAL*8 bhspin_1_state(15,2),bhspin_2_state(15,2)
      REAL*8 dtp_state(15)
      INTEGER binstate_state(15,2)
      COMMON /checkstate_params/ mass_1_state,mass_2_state,
     &           kstar_1_state,kstar_2_state,
     &           sep_state,porb_state,ecc_state,
     &           RROL_1_state,RROL_2_state,evol_type_state,
     &           aj_1_state,aj_2_state,tms_1_state,tms_2_state,
     &           massc_1_state,massc_2_state,rad_1_state,rad_2_state,
     &           mass0_1_state,mass0_2_state,lum_1_state,lum_2_state,
     &           radc_1_state,radc_2_state,menv_1_state,
     &           menv_2_state,renv_1_state,renv_2_state,
     &           omega_spin_1_state,omega_spin_2_state,
     &           B0_1_state,B0_2_state,bacc_1_state,bacc_2_state,
     &           tacc_1_state,tacc_2_state,epoch_1_state,epoch_2_state,
     &           bhspin_1_state,bhspin_2_state,
     &           dtp_state,binstate_state
      REAL*8 checkstate_array(15,123)
      COMMON /checkstate_array/ checkstate_array
