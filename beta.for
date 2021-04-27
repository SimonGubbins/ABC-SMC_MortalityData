      FUNCTION BETA_RAND(NU1,NU2)
C
C--------------------------------------------------------------------

C Function to return a random variable drawn from a 
C beta(NU1,NU2) distribution

C If X1 and X2 are independent random variables, gamma distributed
C with parameters (NU1,1) and (NU2,1) respectively, then

C Y = X1/(X1 + X2)

C is a random variable from the beta(Nu1,Nu2) distribution

C Based on ZQLBET1 from http://www.homepages.ucl.ac.uk/~ucakarc/work/randgen.html
C Uses function gamdev.for which is based on ZBQLGAM from same source

C Note: reparameterisation of gamdev from GAMDEV(G,H) to GAMDEV(shape,mean) where
C G=shape and H=shape/mean

C Tom Sumner - 09/02/2012
C
C--------------------------------------------------------------------
C
      IMPLICIT NONE

      REAL*8 BETA_RAND,NU1,NU2,GAMDEV,X1,X2,dummy

	dummy=1.0

	X1=GAMDEV(NU1,NU1)
	X2=GAMDEV(NU2,NU2)
	BETA_RAND=X1/(X1+X2)

      RETURN
      END      
            