      FUNCTION HOWMANY(N,P)
C
C Function for stochastic simulations. It determines the number of
C events that occur, X, for a population of size N when the individual
C probability for an event is P. This follows a Binomial(N,P)
C distribution, but for certain values of N and P, it can be
C approximated by a Poisson or Normal distribution.
C
C
      IMPLICIT NONE
C
      INTEGER*4 N,HOWMANY
      REAL*8 P
      REAL*8 Z
      REAL*8 MEAN,VAR,NQ
C
C
C Calculate the mean and variance for the Binomial(N,P) distribution
C
      MEAN=DFLOAT(N)*P
      VAR=DFLOAT(N)*P*(1.0-P)
      NQ=DFLOAT(N)*(1.0-P)
C
C Determine if either the Normal or Poisson approximation holds (see
C Evans, Hastings and Peacock (Statistical Distributions)). If not use
C the Binomial distribution
C
      IF ((P.LT.0.1).AND.(MEAN.LT.10.0)) THEN
C
C Use the Poisson approximation
C
        CALL POISSON(HOWMANY,MEAN)

      ELSE IF ((VAR.GT.25.0).OR.
     +         ((VAR.GT.5.0).AND.((P.GE.0.1).AND.(P.LE.0.9))).OR.
     +         ((MEAN.GT.10.0).AND.(NQ.GT.10.0))) THEN

C Use the Normal approximation
C
        CALL NORMAL(Z)
        HOWMANY=INT(MEAN+SQRT(VAR)*Z)
      ELSE
C
C Use the Binomial distribution
C
        CALL BINOMIAL(HOWMANY,N,P)
      ENDIF
C
C Ensure that the normal or poisson approximations have not under- or over-flowed 
C (i.e. we don't get more transitions than there are individuals in the population)
	IF (HOWMANY.LT.0) THEN
          HOWMANY=0
	ELSE IF (HOWMANY.GT.N) THEN
		HOWMANY=N
	ENDIF
C
C
      RETURN
      END            