***
        SUBROUTINE WRITEBPP(jp,tphys,evolve_type,
     &                      mass,kstar,sep,
     &                      tb,ecc,rrl1,rrl2,bkick)
        IMPLICIT NONE
        INCLUDE 'const_bse.h'
*
* Concatenate Strings.
*
*     Author : Scott Coughlin
*     Date :   12th March 2019
*
        REAL*8 bkick(20),mass(2)
        REAL*8 evolve_type,sep,tb,ecc,tphys,rrl1,rrl2
        INTEGER jp,jj
        INTEGER kstar(2)

        jp = MIN(80,jp + 1)
        bpp(jp,1) = tphys
        bpp(jp,2) = mass(1)
        bpp(jp,3) = mass(2)
        bpp(jp,4) = float(kstar(1))
        bpp(jp,5) = float(kstar(2))
        bpp(jp,6) = sep
        bpp(jp,7) = tb
        bpp(jp,8) = ecc
        bpp(jp,9) = rrl1
        bpp(jp,10) = rrl2
        bpp(jp,11) = evolve_type
        bpp(jp,12) = bkick(15)
        bpp(jp,13) = bkick(16)
        if(bkick(1).gt.0.d0.and.bkick(5).le.0.d0)then
            bpp(jp,14) = bkick(13)
            bpp(jp,15) = bkick(18)
        elseif(bkick(1).gt.0.d0.and.bkick(5).gt.0.d0)then
            bpp(jp,14) = bkick(14)
            bpp(jp,15) = bkick(19)
        endif
        DO jj = 13,20
            bkick(jj) = 0.0
        ENDDO

        END
