      FUNCTION GAMDEV(G,H)

C     Returns a random number with a gamma distribution with mean
C     G/H and variance G/(H^2). (ie. shape parameter G & scale
C     parameter H)

C	Modified from routine ZBQLGAM taken from http://www.homepages.ucl.ac.uk/~ucakarc/work/software/randgen.f

C	Tom Sumner 02/11/10

      REAL*8 C,D,R,G,H,A,z1,z2,B1,B2,M
      REAL*8 U1,U2,U,V,TEST,X
      REAL*8 c1,c2,c3,c4,c5,w

	REAL*8 RAN3,GAMDEV
	INTEGER*4 IDUM
C
      COMMON IDUM

      GAMDEV = 0.0

      IF ( (G.LE.0.0).OR.(H.LT.0.0) ) THEN
		WRITE(*,1)
		RETURN
      ENDIF

      IF (G.LT.1.0) THEN
889	u=RAN3(IDUM)
	v=RAN3(IDUM)
	if (u.gt.exp(1.0)/(g+exp(1.0))) goto 891
	GAMDEV=((g+exp(1.0))*u/exp(1.0))**(1.0/g)
	if (v.gt.exp(-GAMDEV)) then
	goto 889
	else
	goto 892
	endif
	WRITE(*,*) 'gamdev1'
891	GAMDEV=-log((g+exp(1.0))*(1.0-u)/(g*exp(1.0)))
	if (v.gt.GAMDEV**(g-1.0)) goto 889
892	GAMDEV=GAMDEV/h
	RETURN
      ELSEIF (G.LT.2.0) THEN
	M = 0.0
	elseif (g.gt.10.0) then
	c1=g-1.0
	c2=(g-1.0/(6.0*g))/c1
	c3=2.0/c1
	c4=c3+2.0
	c5=1.0/sqrt(g)
777	u=RAN3(IDUM)
	v=RAN3(IDUM)
	if (g.gt.2.50) then
	u=v+c5*(1.0-1.860*u)
       endif 
       if (u.le.0.0.or.u.ge.1.0) goto 777 
       w=c2*v/u 
       if (c3*u+w+1.0/w.le.c4) goto 778 
       if (c3*log(u)-log(w)+w.ge.1.0) goto 777
778    GAMDEV=c1*w/h 
       return
      ELSE
       M = -(G-2.0) 
      ENDIF
      R = 0.50
      a = ((g-1.0)/exp(1.0))**((g-1.0)/(r+1.0))
      C = (R*(M+G)+1.0)/(2.0*R)
      D = M*(R+1.0)/R
      z1 = C-sqrt(C*C-D)
*
*     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
*     is not exactly zero - this needs trapping if negative.
*
      IF ((Z1-M.LT.0.0).AND.(Z1-M.GT.-1.0-12)) Z1 = M
      z2 = C+sqrt(C*C-D)
      B1=(z1*(z1-M)**(R*(G-1.0)/(R+1.0)))*EXP(-R*(z1-M)/(R+1.0))
      B2=(z2*(z2-M)**(R*(G-1.0)/(R+1.0)))*EXP(-R*(z2-M)/(R+1.0))
50    U1=RAN3(IDUM)
      U2=RAN3(IDUM)
      U=A*U1
      V=B1+(B2-B1)*U2
      X=V/(U**R)
      IF (X.LE.M) GOTO 50
      TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0))
      IF (U.LE.TEST) THEN
       GAMDEV = (X-M)/H
      ELSE
       GOTO 50
      ENDIF
 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' GAMDEV',/5X, '(both parameters must be positive)',/)

      END
