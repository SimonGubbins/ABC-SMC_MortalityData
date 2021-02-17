	PROGRAM ImplementABC_SMC
C---------------------------------------------------------------------
C ABC_SMC implements an approximate Bayesian sequential Monte Carlo
C sampler to estimate parameters for from mortality data for a high-
C mortality disease
C
C--------------------------------------------------------------------
C DECLARE VARIABLES
C--------------------------------------------------------------------
	IMPLICIT NONE	
C
C DATA---------------------------------------------------------------
C Herd to fit
      INTEGER*4 HFit
C
C Population size
      INTEGER*4 NAnim
C
C Mortality data to fit to
      INTEGER*4 NObs
      INTEGER*4, ALLOCATABLE :: TObs(:)
      INTEGER*4, ALLOCATABLE :: ObsDeadAnim(:)
      INTEGER*4, ALLOCATABLE :: SimDeadAnim(:)
      INTEGER*4 DeadOnFarm
C
C Cumulative mortality
      INTEGER*4 CumMort
C
C Days of first and last observations in the mortality data
      INTEGER*4 T1,T2
C
C Day of confirmation for the population (i.e. day in the mortality
C data when infection confirmed, T1 <= TConf <= T2)
      INTEGER*4 TConf
C
C ABC STUFF----------------------------------------------------------
C Number of parameters
      INTEGER*4, PARAMETER :: nPar = 7
C Number of particles and rounds
	INTEGER*4 nP,nT
      INTEGER*4 T,P,PSel
C Acceptance tolerance
      REAL*8 eps,eps0
C Parameters
      REAL*8, ALLOCATABLE :: theta(:,:,:)
      REAL*8 thetaD(nPar)
C Particle weights
      REAL*8, ALLOCATABLE :: W(:,:)
      REAL*8 sumW
C Perturbation kernel
      REAL*8 Q,s(nPar),UNIFPDF
C Prior stuff
      REAL*8 prior
      INTEGER*4 priorID
      INTEGER*4 priorFlag(nPar-1)
C Fit metric
      REAL*8, ALLOCATABLE :: Dist(:)
C
C MODEL STUFF--------------------------------------------------------
C Time at start of outbreak
      INTEGER*4 T0
C
C Epidemiological parameters:
C kE, muE - shape and mean for latent period
C kI, muI - shape and mean for infectious period
C beta - transmission rate
	INTEGER*4 kE,kI
	REAL*8 muE,muI
      REAL*8 beta
C
C Baseline mortality rate
	REAL*8 rM
C
C Initial number of infected animals
      INTEGER*4 I0
C
C MISCELLANY---------------------------------------------------------
	INTEGER*4 I,J
	INTEGER*4 IDUM
	REAL*8 RAN3,PROB,GAMDEV
C
      COMMON IDUM
C
C--------------------------------------------------------------------
C PREPARE THE INPUTS FOR THE SMC ALGORITHM
C--------------------------------------------------------------------
C Initialise the SMC algorithm
      OPEN(11,FILE='.\SMCInitFile.txt',SHARE='DENYNONE')
	WRITE(*,*) 'Initializing SMC algorithm'
	WRITE(*,*) ' '
C
C Set the random number generator seed and initialise the RNG
      READ(11,*) IDUM
	PROB=RAN3(-IDUM)
C
C Set the population size
      READ(11,*) NAnim
C
C Set the number of observations (i.e. rows in MortalityData.txt)
      READ(11,*) NObs
C
C Set the day of confirmation
      READ(11,*) TConf
C
C Set the initial number of infected animals
      READ(11,*) I0
C
C Set the priors to use
      READ(11,*) priorID
C
C Set the number of particles and number of rounds
      READ(11,*) nP
      READ(11,*) nT
C
C Initial acceptance tolerance
      READ(11,*) eps0
C
      CLOSE(11)
C
C      
C Specify the appropriate prior flags
C Informative all      
      IF (priorID.EQ.1) THEN
          priorFlag(1)=1
          priorFlag(2)=1
          priorFlag(3)=1
          priorFlag(4)=1
          priorFlag(5)=1
          priorFlag(6)=1
C
C Informative all, except transmission
      ELSE IF (priorID.EQ.2) THEN
          priorFlag(1)=1
          priorFlag(2)=1
          priorFlag(3)=1
          priorFlag(4)=1
          priorFlag(5)=0
          priorFlag(6)=1
C
C Informative all, except transmission and mortality
      ELSE IF (priorID.EQ.3) THEN
          priorFlag(1)=1
          priorFlag(2)=1
          priorFlag(3)=1
          priorFlag(4)=1
          priorFlag(5)=0
          priorFlag(6)=0
C
C Informative transmission, otherwise non-informative
      ELSE IF (priorID.EQ.4) THEN
          priorFlag(1)=0
          priorFlag(2)=0
          priorFlag(3)=0
          priorFlag(4)=0
          priorFlag(5)=1
          priorFlag(6)=0
C
C Informative mortality, otherwise non-informative
      ELSE IF (priorID.EQ.5) THEN
          priorFlag(1)=0
          priorFlag(2)=0
          priorFlag(3)=0
          priorFlag(4)=0
          priorFlag(5)=0
          priorFlag(6)=1
C
C Non-informative all
      ELSE IF (priorID.EQ.6) THEN
          priorFlag(1)=0
          priorFlag(2)=0
          priorFlag(3)=0
          priorFlag(4)=0
          priorFlag(5)=0
          priorFlag(6)=0
      END IF
C
C--------------------------------------------------------------------
C LOAD DATA
C--------------------------------------------------------------------
C Create the arrays for the data
      ALLOCATE(TObs(NObs),ObsDeadAnim(NObs))
C
C Load the mortality data
      OPEN(11,FILE='.\MortalityData.txt',SHARE='DENYNONE')
      DO J=1,NObs
	    READ(11,*) TObs(J),ObsDeadAnim(J)
	END DO
      CLOSE(11)
C
C Calculate the cumulative mortality
      CumMort=SUM(ObsDeadAnim)
C
C Set the times of the first and last observation
      T1=MINVAL(TObs)
      T2=MAXVAL(TObs)
C
C Create arrays
      ALLOCATE(SimDeadAnim(T2-T1+1))
C
C--------------------------------------------------------------------
C INITIALISE THE SMC SCHEME
C--------------------------------------------------------------------
C Create the parameter array, weights array and distance metric
      ALLOCATE(theta(nPar,nP,nT))
      ALLOCATE(W(nP,nT))
      ALLOCATE(Dist(nP))
C
C Initialise the weights and parameter arrays
      DO T=1,nT
          DO P=1,nP
              W(P,T)=0.0
              DO J=1,nPar
                  theta(J,P,T)=0.0
              END DO
          END DO
      END DO
C
C Set the initial acceptance tolerance
      eps=eps0
C
C--------------------------------------------------------------------
C IMPLEMENT ABC SMC
C--------------------------------------------------------------------
C For each round ...
      DO T=1,nT
C
C Output the acceptance tolerance
          OPEN(99,FILE='.\Outputs\tolerance.txt',ACCESS='APPEND')
          WRITE(99,'(I6,2X,F12.5)') T,eps
      	CLOSE(99)
C
C For each particle ...
          DO P=1,nP
C              
      		WRITE(*,'(2(A13,I8))') 'Round: ',T,';  Particle: ',P
C
C--------------------------------------------------------------------
C UPDATE THE PARTICLE
C--------------------------------------------------------------------
C For the first round
10            IF (T.EQ.1) THEN
C
C Sample from the priors (0: non-informative; 1: informative)
C Mean latent period
                  IF (priorFlag(1).EQ.0) THEN
                      thetaD(1)=DBLE(5.0)*RAN3(IDUM)
                  ELSE IF (priorFlag(1).EQ.1) THEN
                      thetaD(1)=GAMDEV(DBLE(2.0),DBLE(1.0))
                  END IF
C
C Latent period shape
                  IF (priorFlag(2).EQ.0) THEN
                      thetaD(2)=DBLE(5.0)*RAN3(IDUM)
                  ELSE IF (priorFlag(2).EQ.1) THEN
                      thetaD(2)=DBLE(5.0)*RAN3(IDUM)
                  END IF
C
C Mean infectious period
                  IF (priorFlag(3).EQ.0) THEN
                      thetaD(3)=DBLE(10.0)*RAN3(IDUM)
                  ELSE IF (priorFlag(3).EQ.1) THEN
                      thetaD(3)=GAMDEV(DBLE(10.0),DBLE(5.0))
                  END IF
C
C Infectious period shape parameter
                  IF (priorFlag(4).EQ.0) THEN
                      thetaD(4)=DBLE(20.0)*RAN3(IDUM)
                  ELSE IF (priorFlag(4).EQ.1) THEN
                      thetaD(4)=DBLE(20.0)*RAN3(IDUM)
                  END IF
C
C Transmission parameter
	            IF (priorFlag(5).EQ.0) THEN
                      thetaD(5)=DBLE(10.0)*RAN3(IDUM)
	            ELSE IF (priorFlag(5).EQ.1) THEN
                      thetaD(5)=GAMDEV(DBLE(1.5),DBLE(1.5))
	            END IF
C
C Mortality rate
                  IF (priorFlag(6).EQ.0) THEN
                      thetaD(6)=DBLE(0.005)*RAN3(IDUM)
                  ELSE IF (priorFlag(6).EQ.1) THEN
                      thetaD(6)=-DBLE(0.0002)*DLOG(RAN3(IDUM))
                  END IF
C
C Time of introduction (always non-informative)
                  thetaD(7)=DFLOAT(T1-30)+
     +                      (DFLOAT(TConf)-
     +                       DFLOAT(T1-30))*RAN3(IDUM)
C
C For subsequent rounds
              ELSE IF (T.GT.1) THEN
C
C Sample according to the particle weights
                  CALL MULTINOMIAL(PSel,W(:,T-1),nP)
C
C Perturb the particle
                  DO J=1,nPar
                      thetaD(J)=theta(J,PSel,T-1)+
     +                          (-s(J)+2.0*s(J)*RAN3(IDUM))
                  END DO
C
              END IF
C
C Ensure the sampled parameter set has a non-zero prior probability
              CALL ComputePrior(prior,priorFlag,thetaD,nPar,T1,TConf)
              IF (prior.EQ.0.0) GOTO 10
C
C--------------------------------------------------------------------
C RUN THE MODEL WITH THE SAMPLED PARAMETERS
C--------------------------------------------------------------------
C Extract the parameters
              muE=thetaD(1)
              kE=CEILING(thetaD(2))
              muI=thetaD(3)
              kI=CEILING(thetaD(4))
              beta=thetaD(5)
              rM=thetaD(6)
              T0=CEILING(thetaD(7))
C
C Simulate the outbreak, computing the simulated mortality data for the herd
              CALL SimulateOutbreak(SimDeadAnim,NAnim,
     +                              I0,T0,T1,T2,
     +                              beta,muE,muI,kE,kI,rM)
C
C Check that mortality doesn't exceed the number of animals and, if it does,
C reject the simulation
              DeadOnFarm=SUM(SimDeadAnim)
              IF (DeadOnFarm.GT.NAnim) THEN
                  WRITE(*,'(A22,I8,1X,I8)') 'Deaths>Animals for herd: ',
     +                                      DeadOnFarm,NAnim
                  GOTO 10
              END IF
C
C--------------------------------------------------------------------
C TEST WHETHER OR NOT TO ACCEPT THE SIMULATION
C--------------------------------------------------------------------
C Compute the distance metric
              Dist(P)=0.0
              DO J=1,T2-T1+1
                  Dist(P)=Dist(P)+((DFLOAT(SimDeadAnim(J))-
     +                              DFLOAT(ObsDeadAnim(J)))**2)/
     +                            CumMort
              END DO
C
C If the simulation is too far from the data, draw another parameter set
              IF (Dist(P).GE.eps) GOTO 10
C
C Otherwise, accept the parameter set
              DO J=1,nPar
                  theta(J,P,T)=thetaD(J)
              END DO
C
C--------------------------------------------------------------------
C UPDATE THE PARTICLE WEIGHT
C--------------------------------------------------------------------
C
C For the first round, particles are weighted equally
              IF (T.EQ.1) THEN
                  W(P,T)=1.0
C
C Otherwise ...
              ELSE IF (T.GT.1) THEN
C
C Compute the prior for the parameter set
                  CALL ComputePrior(prior,priorFlag,theta(:,P,T),nPar,
     +                              T1,TConf)
C
C Compute the denominator for the weights
                  sumW=0.0
                  DO I=1,nP
                      Q=1.0
                      DO J=1,nPar
                          Q=Q*UNIFPDF(theta(J,P,T)-theta(J,I,T-1),
     +                                -s(J),s(J))
                      END DO
                      sumW=sumW+W(I,T-1)*Q
                  END DO
C
C Compute the particle weight
                  W(P,T)=prior/sumW
C
              END IF
C
C--------------------------------------------------------------------
C SAVE THE RESULTS
C--------------------------------------------------------------------
C Output the sampled parameters
		    OPEN(99,FILE='.\Outputs\parsamp.txt',ACCESS='APPEND')
		    WRITE(99,'(2(2X,I8),2(2X,F20.7,2X,I8),2(2X,F20.7),2X,I8)')
     +        T,P,muE,kE,muI,kI,beta,rM,T0
		    CLOSE(99)
C
C Output the simulated mortalities
	        OPEN(88,FILE='.\Outputs\Mortality.txt',ACCESS='APPEND')
              WRITE(88,'(26(2X,I8))') T,P,(SimDeadAnim(J),J=1,T2-T1+1)
              CLOSE(88)
C
          END DO    !particle
C
C--------------------------------------------------------------------
C NORMALISE THE PARTICLE WEIGHTS
C--------------------------------------------------------------------
C Compute the ranges for the perturbation kernel 
          DO J=1,nPar
              s(J)=0.1*(MAXVAL(theta(J,:,T))-MINVAL(theta(J,:,T)))
          END DO
C
C Compute the next tolerance for accepting simulations
          I=IDINT(DBLE(0.5)*DFLOAT(nP))
          CALL RealSort(nP,Dist)
          eps=0.5*(Dist(I)+Dist(I+1))
C
C Compute the normalisation constant
          sumW=0.0
          DO J=1,nP
              sumW=sumW+W(J,T)
          END DO
C
C Update the weights
          DO J=1,nP
              W(J,T)=W(J,T)/sumW
          END DO
C
C Output the weights
          OPEN(99,FILE='.\Outputs\weights.txt',ACCESS='APPEND')
          DO P=1,nP
              WRITE(99,'(2(I6,2X),F12.8)') T,P,W(P,T)
          END DO
      	CLOSE(99)
C
C Output theta
          OPEN(99,FILE='.\Outputs\theta.txt',ACCESS='APPEND')
          DO P=1,nP
              WRITE(99,'(I6,2X,I6,31(2X,F20.7))') T,P,
     +        (theta(J,P,T),J=1,nPar)
          END DO
      	CLOSE(99)
C
      END DO   !round
C
C
	END