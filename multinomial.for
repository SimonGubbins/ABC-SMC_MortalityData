C Tom Sumner 08/06/2010
C
	SUBROUTINE MULTINOMIAL(X,P,K)
C
C Subroutine to return X, a single random variable drawn from a 
C Multinomial distribution
C 
C The outcome of the trial is one of K possible outcomes
C The probability of outcome i is P_i and the sum of P_i over [1,K] is 1
C
C The routine needs a Uniform(0,1) random number generator. This code
C assumes there is a function RAN3 (copied from Numerical Recipes)
C which has been initialised and is somewhere the subroutine can find
C it
C
C Might be quicker if we sorted P in descending order first!
C
	IMPLICIT NONE
C
	INTEGER*4 IDUM,K
	INTEGER*4 X
C	
	REAL*8 U,RAN3
	REAL*8 P(K)
	REAL*8 SumP
C
      COMMON IDUM

	SumP=P(1)
	X=1
	U=RAN3(IDUM)
	DO WHILE ((SumP<U).AND.(X.LT.K))
		X=X+1
		SumP=SumP+P(X)
	END DO	
C
	RETURN
	END