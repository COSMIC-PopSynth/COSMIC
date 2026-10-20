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
*     Date : 15th November 2023, 13th October 2026
*
* How to add new timestep_conditions variables:
*
* 1. Adjust array sizes in checkstate.h:
*    checkstate_array (15, x) -> (15, x + 3 * n_new_vars)
*    checkstate_exclude (15, y) -> (15, y + n_new_vars)
* 2. Adjust ncol below
* 3. Add new variable to subroutine definition and current_state_array
* 4. In evolv2.f add new vars to every call of checkstate
* 5. In checkstate.py add new vars to CHECKSTATE_COLUMNS
*
* How conditions are stored (these arrays are filled in by
* set_checkstates in checkstate.py):
*
* Each row jj (1-15) is one sublist of timestep_conditions, e.g.
* ['1<kstar_1<4', 'dtp=0.01'], and dtp_state(jj) is the dtp to use
* when every condition in that row is met. Unused rows have
* dtp_state(jj) = -1.
*
* Each column of current_state_array (e.g. kstar_1) has three slots in
* checkstate_array: lower bound, flag, upper bound. For column number
* param_index these are at 3*(param_index-1) + 1, + 2 and + 3.
*
* The flag says how to compare against the bounds. It is the sum of
*   1 - lower bound is strict, use .gt. (otherwise .ge.)
*   2 - upper bound is strict, use .lt. (otherwise .le.)
*   4 - value must also not equal checkstate_exclude(jj, param_index)
* Both bounds are always checked, so a condition that only gives one
* side (e.g. kstar_1>1) leaves the other at its default of -10e30 or
* 10e30, which every value passes.
*
* Examples for kstar_1:
*   '1<kstar_1<4'   -> [1, 3, 4]         (1 + 2, both strict)
*   '1<=kstar_1<4'  -> [1, 2, 4]         (only the upper is strict)
*   'kstar_1=4'     -> [4, 0, 4]         (4 <= kstar_1 <= 4)
*   'kstar_1>1'     -> [1, 1, 10e30]
*   'kstar_1!=4'    -> [-10e30, 4, 10e30] and checkstate_exclude = 4
*   no condition    -> [-10e30, 0, 10e30] (column is skipped)
*
* A row passes if every column with a condition passes, so a row with
* only a dtp (no conditions) always passes. If several rows pass then
* the last one sets dtp. If no rows pass then dtp goes back to
* dtp_original (the dtp the binary was evolved with).
*
      IMPLICIT NONE
      INCLUDE 'checkstate.h'
      INCLUDE 'const_bse.h'
*  number of columns that can be used as timestep conditions
      INTEGER ncol
      PARAMETER (ncol = 43)
      INTEGER jj,ii,param_index,binstate,iflag
      INTEGER kstar1,kstar2
      LOGICAL pass_condition,pass_condition_any
      LOGICAL isave,iplot
      REAL*8 current_state_array(ncol)
      REAL*8 lo,hi,val,dtp_prev
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

* tsave is the next time at which a line will be written to the bcm
* array, it should never be bigger than tphysf
      IF(tsave.ge.tphysf)THEN
          tsave = tphysf
      ENDIF

* Remember the dtp used up to now, so that once we have worked out the
* new dtp we can tell whether it has changed (see the end of the routine)
      dtp_prev = dtp

* We track to see if any row passes, otherwise dtp will be set to dtp_original
      pass_condition_any = .false.

* Loop over each row of conditions
      DO jj = 1,15
* First check if there is even a dtp requested, if not
* this implies that this row is unused so we skip to the next one
          IF(dtp_state(jj).eq.-1.d0)THEN
              goto 74
          ENDIF
* Assume the row passes until one of its conditions fails. This means
* that a row with only a dtp (and no conditions) always passes
          pass_condition = .true.

* Loop over each column, checking its condition (if it has one)
          DO param_index = 1,ncol
* Index of the lower bound slot for this column (the flag and upper
* bound are in the next two slots)
              ii = 3*(param_index-1) + 1
              lo = checkstate_array(jj,ii)
* The flag is stored as a real (all slots share one array) so round it
* back to an integer so that we can test its bits with IAND
              iflag = NINT(checkstate_array(jj,ii+1))
              hi = checkstate_array(jj,ii+2)

* If the slots still have their default values then this column has no
* condition in this row, so skip it. Without this skip, the bounds
* checks below would still pass for these defaults, but any column with
* a NaN value would make every row fail
              IF(iflag.eq.0.and.lo.eq.-10.d30.and.hi.eq.10.d30)THEN
                  goto 69
              ENDIF
              val = current_state_array(param_index)

* Check the lower bound, strict (val > lo) if bit 1 is set, otherwise
* val >= lo. The tests are written as .not.(passes) rather than
* (fails) because any comparison with a NaN is false, so this way
* a NaN value fails the condition rather than passing it
              IF(IAND(iflag,1).ne.0)THEN
                  IF(.not.(val.gt.lo)) pass_condition = .false.
              ELSE
                  IF(.not.(val.ge.lo)) pass_condition = .false.
              ENDIF

* Check the upper bound in the same way, strict (val < hi) if bit 2 is
* set, otherwise val <= hi
              IF(IAND(iflag,2).ne.0)THEN
                  IF(.not.(val.lt.hi)) pass_condition = .false.
              ELSE
                  IF(.not.(val.le.hi)) pass_condition = .false.
              ENDIF

* If bit 4 is set then this column also has a != condition, so the
* row fails if the value equals the excluded one
              IF(IAND(iflag,4).ne.0)THEN
                  IF(val.eq.checkstate_exclude(jj,param_index))THEN
                      pass_condition = .false.
                  ENDIF
              ENDIF

* As soon as one condition is false the row has failed, so we do not
* check the rest of its columns
              IF(.not.pass_condition)THEN
                  goto 79
              ENDIF
 69           CONTINUE
          ENDDO
 79       CONTINUE

* If the row passed then use its dtp. Later rows are checked after this
* one and overwrite dtp if they pass too, so the last passing row wins
          IF(pass_condition)THEN
              dtp = dtp_state(jj)
              pass_condition_any = .true.
          ENDIF
 74   CONTINUE
      ENDDO

* If no rows passed then return to the original dtp
      IF(pass_condition_any.neqv..true.)THEN
          dtp = dtp_original
      ENDIF

* If dtp has changed then tsave (the next output time) may still be set
* by the old dtp, since evolv2.f only updates it after writing a line
* (tsave = tsave + dtp). Two cases need fixing:
*   1. switching to a finer dtp, e.g. from 50 to 0.01 Myr, tsave could
*      still be up to 50 Myr away, missing short phases entirely
*   2. switching away from dtp = 0, tsave was set to tphysf (below) so
*      no more lines would be written until the end of the evolution
* So we bring tsave forward to at most one new dtp from now. Taking the
* MIN means that when switching to a coarser dtp we keep the save time
* that was already scheduled. dtp == 0 is handled separately below
      IF(dtp.ne.dtp_prev.and.dtp.gt.0.d0)THEN
          tsave = MIN(tsave,tphys + dtp)
      ENDIF

* Set how evolv2.f writes to the bcm array: isave means write a line
* once tphys reaches tsave, while iplot means write a line every
* timestep. For dtp == 0 we write every timestep, so tsave is moved to
* tphysf to stop it limiting the timestep size
      isave = .true.
      iplot = .false.
      IF(dtp.eq.0.d0)THEN
         iplot = .true.
         isave = .false.
         tsave = tphysf
      ENDIF
      END
