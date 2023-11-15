***
      SUBROUTINE checkstate(dtp,dtp_original,tsave,tphys,tphysf,
     &                      iplot,isave,binstate,evolve_type,
     &                      mass1,mass2,kstar1,kstar2,sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj1,aj2,tms1,tms2,
     &                      massc1,massc2,rad1,rad2,
     &                      mass0_1,mass0_2,lumin1,lumin2,
     &                      radc1,radc2,menv1,menv2,renv1,renv2,
     &                      ospin1,ospin2,b_0_1,b_0_2,bacc1,bacc2,
     &                      tacc1,tacc2,epoch1,epoch2,bhspin1,bhspin2,
     &                      teff1,teff2)
*
* Check timestep conditions for bcm array
*
*     Author : Scott Coughlin
*     Date :   7th April 2020
*
*     Edited : Tom Wagg
*     Date : 15th November 2023
*
* How to add new timestep_conditions variables:
*
* 1. Adjust array size of checkstate_array in checkstate.h (15, x) -> (15, x + 3 * n_new_vars)
* 2. Adjust current_state_array size below
* 3. Add new variable to subroutine definition and current_state_array
* 4. Change DO loop conditions to (2, x, 3) -> (2, x + 3 * n_new_vars, 3)
* 5. In evolv2.f add new vars to every call of checkstate
* 6. In checkstate.py add new vars to CHECKSTATE_COLUMNS
*
      IMPLICIT NONE
      INCLUDE 'checkstate.h'
      INCLUDE 'const_bse.h'
      INTEGER jj,ii,param_index,binstate
      INTEGER kstar1,kstar2
      LOGICAL pass_condition,pass_condition_any
      LOGICAL isave,iplot
*  current_state_array length set by number of columns that can be used as timestep conditions
      REAL*8 current_state_array(43)
      REAL*8 dtp,dtp_original,tsave,tphys,tphysf,mass1,mass2
      REAL*8 evolve_type,sep,tb,ecc,rrl1,rrl2
      REAL*8 aj1,aj2,tms1,tms2,massc1,massc2,rad1,rad2
      REAL*8 mass0_1,mass0_2,lumin1,lumin2,radc1,radc2
      REAL*8 menv1,menv2,renv1,renv2,ospin1,ospin2
      REAL*8 b_0_1,b_0_2,bacc1,bacc2,tacc1,tacc2,epoch1,epoch2
      REAL*8 bhspin1,bhspin2,teff1,teff2
      current_state_array(1) = binstate
      current_state_array(2) = evolve_type
      current_state_array(3) = mass1
      current_state_array(4) = mass2
      current_state_array(5) = kstar1
      current_state_array(6) = kstar2
      current_state_array(7) = sep
      current_state_array(8) = tb
      current_state_array(9) = ecc
      current_state_array(10) = rrl1
      current_state_array(11) = rrl2
      current_state_array(12) = aj1
      current_state_array(13) = aj2
      current_state_array(14) = tms1
      current_state_array(15) = tms2
      current_state_array(16) = massc1
      current_state_array(17) = massc2
      current_state_array(18) = rad1
      current_state_array(19) = rad2
      current_state_array(20) = mass0_1
      current_state_array(21) = mass0_2
      current_state_array(22) = lumin1
      current_state_array(23) = lumin2
      current_state_array(24) = radc1
      current_state_array(25) = radc2
      current_state_array(26) = menv1
      current_state_array(27) = menv2
      current_state_array(28) = renv1
      current_state_array(29) = renv2
      current_state_array(30) = ospin1
      current_state_array(31) = ospin2
      current_state_array(32) = b_0_1
      current_state_array(33) = b_0_2
      current_state_array(34) = bacc1
      current_state_array(35) = bacc2
      current_state_array(36) = tacc1
      current_state_array(37) = tacc2
      current_state_array(38) = epoch1
      current_state_array(39) = epoch2
      current_state_array(40) = bhspin1
      current_state_array(41) = bhspin2
      current_state_array(42) = teff1
      current_state_array(43) = teff2

* tsave should never be bigger than tphysf
      IF(tsave.ge.tphysf)THEN
          tsave = tphysf
      ENDIF

* We track to see if any conditional pass, otherwise dtp will be set to
* dtp_original
      pass_condition_any = .false.

      DO jj = 1,15
* start by assuming we do not pass the current conditional setting of dtp
          pass_condition = .false.
* First check if there is even a dtp requested, if not
* this implies that a conditional has not been set so we should continue
          IF(dtp_state(jj).eq.-1.d0)THEN
              goto 74
          ENDIF
* if we do have the conditional then we need to know the sign of the conditional
* i.e. EQUAL (0), GT (1), GE (2), LT (3), LE (4) and then piece wise together the whole statement
* Moreover we need to keep track of which parameter in the current_state_array we are checking
          param_index = 0
* Loop stop value here is (number of condition columns - 1) * 3
          DO ii = 2, 126, 3
              param_index = param_index + 1
* This part of the checkstate array had no conditional set (we know because the default has not changed)
              IF(checkstate_array(jj,ii-1).eq.-10E30.and.
     &               checkstate_array(jj,ii+1).eq.10E30)THEN
                  goto 69
              ENDIF

* Now we check for what the conditional was and apply it. As soon as the
* conditional is false we do not check anymore
              IF(checkstate_array(jj,ii).eq.0)THEN

                  IF(checkstate_array(jj,ii-1).eq.
     &              current_state_array(param_index)
     &              .and.
     &              checkstate_array(jj,ii+1).eq.
     &              current_state_array(param_index))THEN
                      pass_condition = .true.
                  ELSE
                      pass_condition = .false.
                      goto 79
                  ENDIF

              ELSEIF(checkstate_array(jj,ii).eq.1)THEN

                  IF(current_state_array(param_index).gt.
     &               checkstate_array(jj,ii-1))THEN
                      pass_condition = .true.
                  ELSE
                      pass_condition = .false.
                      goto 79
                  ENDIF
              
              ELSEIF(checkstate_array(jj,ii).eq.2)THEN
                  IF(current_state_array(param_index).ge.
     &               checkstate_array(jj,ii-1))THEN
                      pass_condition = .true.
                  ELSE
                      pass_condition = .false.
                      goto 79
                  ENDIF

              ELSEIF(checkstate_array(jj,ii).eq.3)THEN

                  IF(current_state_array(param_index).lt.
     &               checkstate_array(jj,ii+1))THEN
                      pass_condition = .true.
                  ELSE
                      pass_condition = .false.
                      goto 79
                  ENDIF

              ELSEIF(checkstate_array(jj,ii).eq.4)THEN
                  IF(current_state_array(param_index).le.
     &               checkstate_array(jj,ii+1))THEN
                      pass_condition = .true.
                  ELSE
                      pass_condition = .false.
                      goto 79
                  ENDIF
              ENDIF
 69       CONTINUE
          ENDDO
 79       CONTINUE
* finally we see if we satisfied the conditional and set dtp
          IF(pass_condition)THEN
              IF(dtp_state(jj).ne.tphysf.and.tsave.ge.tphysf)THEN
                  tsave = tphys
              ENDIF
              dtp = dtp_state(jj)
              pass_condition_any = .true.
          ENDIF
 74   CONTINUE
      ENDDO
      IF(pass_condition_any.neqv..true.)THEN
* if no conditions are met then return to orginal dtp
          dtp = dtp_original
      ENDIF
      isave = .true.
      iplot = .false.
      IF(dtp.eq.0.d0)THEN
         iplot = .true.
         isave = .false.
         tsave = tphysf
      ENDIF
      END
