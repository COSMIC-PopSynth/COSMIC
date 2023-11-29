*
* checkstate.h
*
      REAL*8 dtp_state(15)
      COMMON /checkstate_params/ dtp_state
      REAL*8 checkstate_array(15,129)               ! (15, 3 * n_col)
      COMMON /checkstate_array/ checkstate_array
      INTEGER check_dtp
      COMMON /check_dtp/ check_dtp
