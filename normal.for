      SUBROUTINE NORMAL(Z)
C
C Subroutine to return Z a random variable drawn from a Normal(0,1)
C distribution
C
C For an outline of the Box-Muller method see Ross (A First Course in
C Probabilty, pg 395-6). This produces two normal variates, so output
C one value and store the other for the next call to the routine.
C
C The routine needs a Uniform(0,1) random number generator. This code
C assumes there is a function RAN3 (copied from Numerical Recipes)
C which has been initialised and is somewhere the subroutine can find
C it
C
C
      IMPLICIT NONE
C
      INTEGER*4 IDUM,IPREV
      REAL*8 RAN3,Z,ZSPARE,V1,V2,S
C
      COMMON IDUM
C
10    V1=2.0*RAN3(IDUM)-1.0
      V2=2.0*RAN3(IDUM)-1.0
C
      S=V1**2+V2**2
      IF ((S.GT.1).OR.(S.EQ.0)) GOTO 10
C
      Z=V1*SQRT(-2.0*LOG(S)/S)
C
      RETURN
      END      
            