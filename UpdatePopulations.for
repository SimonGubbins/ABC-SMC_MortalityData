	SUBROUTINE UpdatePopulations(DT,S,E,I,kE,kI,
     +						     beta,muE,muI,rM,NewDeadAnim)
C--------------------------------------------------------------------
C Subroutine to update the populations from one time step to the next
C for the stochastic epidemic model for a directly transmitted high-
C mortality disease (i.e. one in which all animals can be assumed to die)
C
C--------------------------------------------------------------------
C DECLARE VARIABLES
C--------------------------------------------------------------------
	IMPLICIT NONE
C
C Time step
	REAL*8 DT
C
C Susceptible (S), latent (E) (in each kE stage), infectious (I) (in each kI stage)
C animals
	INTEGER*4 kE,kI
	INTEGER*4 S,E(kE),I(kI)
      INTEGER*4 NAnim
C
C Number of transitions (infection, death, advance through infection
C stages and recruitment)
	INTEGER*4 NInf,NAdvE(kE),NAdvI(kI)
	INTEGER*4 NDeadS,NDeadE(kE),NDeadI(kI)
      INTEGER*4 NewDeadAnim
C
C Force of infection
	REAL*8 Lambda
C
C Parameters
	REAL*8 muE,muI
      REAL*8 rE,rI,rM
	REAL*8 beta
C
C Miscellany
	INTEGER*4 J,HOWMANY
C
C Compute rates for progress through each infection stage
      rI=1/muI
      rE=1/muE
C
C--------------------------------------------------------------------
C COMPUTE FORCE OF INFECTION
C--------------------------------------------------------------------
	NAnim=S+SUM(E)+SUM(I)
	Lambda=beta*(DFLOAT(SUM(I))/DFLOAT(NAnim))
C
C--------------------------------------------------------------------
C COMPUTE THE NUMBER OF TRANSITIONS OF EACH TYPE
C--------------------------------------------------------------------
C
C Susceptible animals
C
C Infection
	NInf=0
	IF ((S.GT.0).AND.(Lambda.GT.0)) THEN
	    NInf=HOWMANY(S,Lambda*DT)
      ENDIF
C
C Natural mortality
      NDeadS=0
	IF ((S.GT.0).AND.(rM.GT.0)) THEN
	    NDeadS=HOWMANY(S,rM*DT)
      ENDIF
C
C Latent animals
      DO J=1,kE
C
C Advance through latent stages
	    NAdvE(J)=0
	    IF ((E(J).GT.0).AND.(rE.GT.0.0)) THEN
	        NAdvE(J)=HOWMANY(E(J),DFLOAT(kE)*rE*DT)
          END IF
C
C Natural mortality
          NDeadE(J)=0
	    IF ((E(J).GT.0).AND.(rM.GT.0.0)) THEN
	        NDeadE(J)=HOWMANY(E(J),rM*DT)
          END IF
C
      END DO
C
C Infectious animals
	DO J=1,kI
C
C Advance through infectious stages
	    NAdvI(J)=0
	    IF ((I(J).GT.0).AND.(rI.GT.0.0)) THEN
	        NAdvI(J)=HOWMANY(I(J),DFLOAT(kI)*rI*DT)
          END IF
C
C Natural mortality
          NDeadI(J)=0
	    IF ((I(J).GT.0).AND.(rM.GT.0.0)) THEN
	        NDeadI(J)=HOWMANY(I(J),rM*DT)
          END IF
C
      END DO
C
C--------------------------------------------------------------------
C UPDATE POPULATIONS (ensuring they are always non-negative)
C--------------------------------------------------------------------
C Susceptible animals
	S=MAX(0,S-NInf-NDeadS)
C
C Latent animals (stage 1)
      E(1)=MAX(0,E(1)+NInf-NAdvE(1)-NDeadE(1))
C
C Latent animals (stages 2 to kE)
	DO J=2,kE
	    E(J)=MAX(0,E(J)+NAdvE(J-1)-NAdvE(J)-NDeadE(J))
      END DO
C
C Infectious animals (stage 1)
	I(1)=MAX(0,I(1)+NAdvE(kE)-NAdvI(1)-NDeadI(1))
C
C Infectious animals (stages 2 to kI)
	DO J=2,kI
	    I(J)=MAX(0,I(J)+NAdvI(J-1)-NAdvI(J)-NDeadI(J))
      END DO
C
C Compute the number of newly dead animals (infected animals are
C assumed to die at the end of the infectious period)
      NewDeadAnim=NDeadS+SUM(NDeadE)+SUM(NDeadI)+NAdvI(kI)
C
C
	RETURN
	END