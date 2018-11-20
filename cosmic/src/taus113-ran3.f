C ====================================================================
C     This function is a drop-in replacement for the ran3 random
C     number generator from Numerical Recipes, Press et al.
C     It is based on E'cuyer's combined Tausworth generator with
C     period 2**113.
C
C ====================================================================
      real function ran3(IDUM)
        implicit none
        include "tausworth.h"
        integer(kind=int64) state(4), rnumber, taus113_gen_int
        integer i, IDUM, first, tid
c$      integer OMP_GET_THREAD_NUM
        COMMON/Taus113State/ state, first
        SAVE /Taus113State/
!$OMP THREADPRIVATE(/Taus113State/)
        EXTERNAL taus113_gen_int

        if(IDUM.lt.0) then
          call taus113_seeding(state, -IDUM)
C         The "warm-up" that has been left out in the seeding routine:
          print*, "Combined Tausworth 113 Generator warming up."
          do i=1,10
            call taus113_next_state(state)
          end do
          IDUM= -IDUM
        endif

        call taus113_next_state(state)
        rnumber= taus113_gen_int(state)

        if (first.eq.1) then
          tid= 1
c$        tid= OMP_GET_THREAD_NUM()
          first=0
          print("('taus-ran3: TID',I2,' The first random number is ',
     &    G16.8)"), tid, real(dble(rnumber)/dble(umax))
        endif

        ran3= real(dble(rnumber)/dble(umax))
      end function

      subroutine taus113_ran3_tester()
        integer i, failcount
        logical passed
        real rnumber, ran3
        external ran3

        passed=.true.
        failcount=0
        rnumber=ran3(-123456)
        do i=1,1000000000
          if (rnumber>=1.e0) then
            passed=.false.
            failcount= failcount + 1
            print*, "FOUL!!", rnumber-1.0e0, rnumber==1.e0, failcount, i
          end if
          rnumber=ran3(123456)
        end do
      end subroutine
