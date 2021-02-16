      SUBROUTINE POISSON(X,LAMBDA)
C
C Subroutine to return X a random variable drawn from a Poisson(LAMBDA)
C distribution (copied from Numerical Recipies)
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
      INTEGER*4 IDUM,X
C
      REAL*8 RAN3,LAMBDA,V,ELAMBDA
      REAL*8 XR,TMP1,TMP2,Y,T
      REAL*8 gammln,pi
      COMMON IDUM
C
      pi=4.0*DATAN(DBLE(1.0))
C
C For small-ish values of LAMBDA simulate a bunch of exponential variables,
C the number of which required to reach a threshold is Poisson
      IF (LAMBDA.LT.DBLE(1000.0)) THEN
C
          V=1.0
	    ELAMBDA=DEXP(-LAMBDA)
          X=-1
10        X=X+1
          V=V*RAN3(IDUM)
          IF (V.GT.ELAMBDA) GOTO 10
C
      ELSE
C
C For larger values of LAMBDA use a rejection method instead
          TMP1=DSQRT(2.0*LAMBDA)
          TMP2=gammln(LAMBDA+1.0)-LAMBDA*DLOG(LAMBDA)
20        Y=DTAN(pi*RAN3(IDUM))
          X=INT(LAMBDA+TMP1*Y)
          IF (X.LT.0) GOTO 20
          XR=DBLE(X)
          T=(XR*DLOG(LAMBDA)-gammln(XR+1.0))+TMP2
          IF (DABS(T).LT.100.0) THEN
              T=0.9*(1.0+(Y*Y))*DEXP(T)
              IF (RAN3(IDUM).GT.T) GOTO 20
          ELSE
              T=DLOG(0.9*(1.0+(Y*Y)))+T
              IF (DLOG(RAN3(IDUM)).GT.T) GOTO 20
          ENDIF
C          
      END IF
C
      RETURN
      END