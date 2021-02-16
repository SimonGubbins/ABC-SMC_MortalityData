      SUBROUTINE BINOMIAL(X,N,P)
C
C Subroutine to return X, a random variable drawn from a Binomial(N,p)
C distribution
C
C For an outline of a routine see Ross (A First Course in Probabilty,
C pg 398)
C
C The routine needs a Uniform(0,1) random number generator. This code
C assumes there is a function RAN3 (copied from Numerical Recipes)
C which has been initialised and is somewhere the subroutine can find
C it
C
C
      IMPLICIT NONE
C
      INTEGER*4 I,IDUM,N,X
C
      REAL*8 P,U,RAN3
C
      COMMON IDUM
C
C
      X=0
      DO I = 1,N
        U=RAN3(IDUM)
        IF (U.LT.P) THEN
          X=X+1
        ENDIF
      END DO
C
      RETURN
      END      
            