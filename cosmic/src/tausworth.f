      BLOCK DATA Taus113_init
        include "tausworth.h"
        DATA c /z'fffffffe', z'fffffff8', z'fffffff0', z'ffffff80'/
        DATA cseed /69069/
      END

      subroutine taus113_seeding(z, seed)
        include "tausworth.h"
        integer(kind=int64) z(4), seed

        print*, umax, bit_size(umax)

        z(1) = iand(umax, seed * cseed)
        if ( z(1) < 2 ) z(1) = z(1) + 2
        z(2) = iand(umax, z(1) * cseed)
        if ( z(2) < 8 ) z(2) = z(2) + 8
        z(3) = iand(umax, z(2) * cseed)
        if ( z(3) < 16 ) z(3) = z(3) + 16
        z(4) = iand(umax, z(3) * cseed)
        if ( z(4) < 128 ) z(4) = z(4) + 128
      end subroutine

      function ishft32(a, shift)
        include "tausworth.h"
        integer(kind=int64) ishft32, a
        integer shift

        ishft32= iand(ishft(a, shift), umax)
      end function

      subroutine taus113_next_state(z)
        include "tausworth.h"
        integer(kind=int64) z(4), znew(4), ishft32
        EXTERNAL ishft32
        integer i

        znew(1)  = ishft32(iand(z(1),c(1)), 18)
        znew(2)  = ishft32(iand(z(2),c(2)), 2)
        znew(3)  = ishft32(iand(z(3),c(3)), 7)
        znew(4)  = ishft32(iand(z(4),c(4)),13)

        znew(1) = ieor(znew(1), ishft(ieor(ishft32(z(1),6), z(1)),-13))
        znew(2) = ieor(znew(2), ishft(ieor(ishft32(z(2),2), z(2)),-27))
        znew(3) = ieor(znew(3), ishft(ieor(ishft32(z(3),13), z(3)),-21))
        znew(4) = ieor(znew(4), ishft(ieor(ishft32(z(4),3), z(4)),-12))

        do i=1,4
          z(i)=znew(i)
        end do
      end subroutine

      function taus113_gen_int(z)
        include "tausworth.h"
        integer(kind=int64) taus113_gen_int, z(4), rnum

        rnum= z(1)
        do i=2,4
          rnum= ieor(rnum, z(i))
        end do

        taus113_gen_int= rnum
      end function

      subroutine test_taus()
        include "int64.h"
        integer(kind=int64) taus113_gen_int
        EXTERNAL taus113_gen_int
        integer(kind=int64) i, z(4), znew(4), seed, maxiter

        seed= 123456
        call taus113_seeding(z, seed)

        maxiter= int(1d9)
        do i=1,maxiter
          call taus113_next_state(z)
        end do

        print*, bit_size(seed), seed, umax
        do i=1,4
          print*, "z(",i,")=", z(i)
        enddo
        print*, taus113_gen_int(z)
      end subroutine

