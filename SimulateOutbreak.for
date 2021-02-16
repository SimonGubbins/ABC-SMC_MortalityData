	SUBROUTINE SimulateOutbreak(NDeadAnim,NAnim,
     +                            I0,T0,T1,TStop,
     +                            beta,muE,muI,kE,kI,rM)
C
C--------------------------------------------------------------------
C SimulateOutbreak simulates within farm dynamics on a farm
C
C Uses a DT=1/100 day timestep
C
C Calls:
C
C 1. UpdatePopulations - updates populations for each time step based
C    on transmission events
C
C--------------------------------------------------------------------
C DECLARE VARIABLES
C--------------------------------------------------------------------
	IMPLICIT NONE
C
C Number of pigs infected initially
	INTEGER*4 I0
C
C Times
	INTEGER*4 T0,T1,TStart,TStop
C
C Simulation time and timestep
      REAL*8 T
	REAL*8, PARAMETER :: DT = 0.01
C
C Flag indicating whether or not infection has been introduced
      INTEGER*4 IFlag
C
C Total, susceptible, infected (in each stage) animals; numbers of
C dead animals
	INTEGER*4 kE,kI
	INTEGER*4 NAnim,S,E(kE),I(kI)
	INTEGER*4 NewDeadAnim
	INTEGER*4 NDeadAnim(TStop-T1+1)
C
C Transmission rate
	REAL*8 beta
C
C Mean duration of latent and infectious periods
	REAL*8 muE,muI
C
C Mortality rate (baseline)
      REAL*8 rM
C
C Miscellany
	INTEGER*4 J
C	
C--------------------------------------------------------------------
C INITIALISE THE SIMULATION
C--------------------------------------------------------------------
C Set times
      TStart=MIN(T0,T1)
	T=DFLOAT(TStart)
      IFlag=0
C
C Set the number of susceptible, exposed and infected animals
	S=NAnim
	DO J=1,kE
	  E(J)=0
      END DO
      E(1)=0
	DO J=1,kI
	  I(J)=0
	END DO
C
C Set the number of dead animals
      DO J=1,TStop-T1+1
          NDeadAnim(J)=0
      END DO
C
C--------------------------------------------------------------------
C SIMULATE WITHIN-FARM DYNAMICS
C--------------------------------------------------------------------
C
	DO WHILE (T.LE.DFLOAT(TStop+1))
C
C Introduce infection
        IF ((FLOOR(T).EQ.T0).AND.(IFlag.EQ.0)) THEN
            S=S-I0
            E(1)=I0
            IFlag=1
        END IF
C
C Update time
	  T=T+DT
C
C Update populations
	  CALL UpdatePopulations(DT,S,E,I,kE,kI,
     +						 beta,muE,muI,rM,NewDeadAnim)
C
C Update the number of dead animals
        J=FLOOR(T)-T1+1
        IF ((J.GT.0).AND.(J.LE.TStop-T1+1)) THEN
            NDeadAnim(J)=NDeadAnim(J)+NewDeadAnim
        END IF
C
	END DO
C
	RETURN
	END