      SUBROUTINE ComputePrior(prior,priorFlag,theta,nPar,T1,TConf)
C
C Subroutine to compute the prior probability for the input parameters
C
C Declare variables
      IMPLICIT NONE
C
C Arguments
      INTEGER*4 nPar,nHerds
      INTEGER*4 T1,TConf
      REAL*8 prior
      REAL*8 theta(nPar)
C
C Flag indicating prior to use for each parameter (0: non-informative, 1: informative)
      INTEGER*4 priorFlag(nPar-1)
C
C Some functions for the PDFs
      REAL*8 UNIFPDF,EXPPDF,GAMPDF
C
C
C Prior for mean latent period
      IF (priorFlag(1).EQ.0) THEN
          prior=UNIFPDF(theta(1),DBLE(0.0),DBLE(5.0))
      ELSE IF (priorFlag(1).EQ.1) THEN
          prior=GAMPDF(theta(1),DBLE(2.0),DBLE(1.0))
      END IF
C
C Prior for latent period shape
      IF (priorFlag(2).EQ.0) THEN
          prior=prior*UNIFPDF(theta(2),DBLE(0.0),DBLE(5.0))
      ELSE IF (priorFlag(2).EQ.1) THEN
          prior=prior*UNIFPDF(theta(2),DBLE(0.0),DBLE(5.0))
      END IF
C
C Prior for mean infectious period parameter
      IF (priorFlag(3).EQ.0) THEN
          prior=prior*UNIFPDF(theta(3),DBLE(0.0),DBLE(10.0))
      ELSE IF (priorFlag(3).EQ.1) THEN
          prior=prior*GAMPDF(theta(3),DBLE(10.0),DBLE(5.0))
      END IF
C
C Prior for infectious period shape parameter
      IF (priorFlag(4).EQ.0) THEN
          prior=prior*UNIFPDF(theta(4),DBLE(0.0),DBLE(20.0))
      ELSE IF (priorFlag(4).EQ.1) THEN
          prior=prior*UNIFPDF(theta(4),DBLE(0.0),DBLE(20.0))
      END IF
C
C Prior for transmission parameter
      IF (priorFlag(5).EQ.0) THEN
          prior=prior*UNIFPDF(theta(5),DBLE(0.0),DBLE(10.0))
      ELSE IF (priorFlag(5).EQ.1) THEN
          prior=prior*GAMPDF(theta(5),DBLE(1.5),DBLE(1.5))
      END IF
C
C Prior for mortality rate
      IF (priorFlag(6).EQ.0) THEN
          prior=prior*UNIFPDF(theta(6),DBLE(0.0),DBLE(0.005))
      ELSE IF (priorFlag(6).EQ.1) THEN
          prior=prior*EXPPDF(theta(6),DBLE(0.0002))
      END IF
C
C Prior for time of introduction
      prior=prior*UNIFPDF(theta(7),DFLOAT(T1-30),DFLOAT(TConf))
C
C      
      RETURN
      END
C
C
C
C------------------------------------------------------
C FUNCTIONS TO COMPUTE PDFs
C------------------------------------------------------
C
C Uniform
      FUNCTION UNIFPDF(X,A,B)
      IMPLICIT NONE
      REAL*8 UNIFPDF,X,A,B
      IF ((X.LT.A).OR.(X.GT.B)) THEN
          UNIFPDF=0.0
      ELSE
          UNIFPDF=1.0/(B-A)
      END IF
      RETURN
      END
C
C Exponential
      FUNCTION EXPPDF(X,mu)
      IMPLICIT NONE
      REAL*8 EXPPDF,X,mu
      IF (X.LT.0.0) THEN
          EXPPDF=0.0
      ELSE
          EXPPDF=(1.0/mu)*DEXP(-X/mu)
      END IF
      RETURN
      END
C
C Gamma
C note: shape=a, mean=a*b
      FUNCTION GAMPDF(X,shape,mean)
      IMPLICIT NONE
      REAL*8 GAMPDF,X,a,b,shape,mean
      REAL*8 G,gammln
      a=shape
      b=mean/shape
      IF (X.LT.0.0) THEN
          GAMPDF=0.0
      ELSE
          G=DEXP(gammln(a))
          GAMPDF=((X/b)**(a-1))*DEXP(-X/b)/(b*G)
      END IF
      RETURN
      END
