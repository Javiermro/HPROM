      SUBROUTINE VECASI(N,V1,V2)
C*****************************************************************************
C
C***  VECTOR ASSIGN:     V2(I) = V1(I)     ->   I=1..N
C
C*****************************************************************************
      IMPLICIT NONE
C
      INTEGER   N
      REAL*8    V1(N),V2(N)
C
      INTEGER   I
C
      DO I=1,N
        V2(I)=V1(I)
      END DO
C
      END
