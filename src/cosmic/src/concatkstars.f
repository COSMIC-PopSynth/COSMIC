***
      SUBROUTINE CONCATKSTARS(int1,int2,int3)
*
* Concatenate Strings.
*
*     Author : Scott Coughlinn
*     Date :   11th December 2018
*
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
      INTEGER int1,int2,int3
      CHARACTER(len=2) str1,str2
      CHARACTER(len=4) str3

      IF(USING_CMC.EQ.0)THEN

          WRITE(str1,'(i2.2)') int1
          WRITE(str2,'(i2.2)') int2

          str3 = trim(adjustl(str1))//trim(adjustl(str2))

          READ(str3,*)int3
      ELSE
          int3 = -1
      ENDIF

      END
