      FUNCTION gammln(xx)
C
C Function to return ln(gamma(xx)), equivalent to ln((xx-1)!) when x is an integer 
C (copied from Numerical Recipes)

C Tom Sumner - 01/11/2010

      IMPLICIT NONE

	REAL*8 gammln,xx
	INTEGER*4 j
	REAL*8 ser,stp,tmp,x,y,cof(6)

	SAVE cof,stp
	DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     +24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     +-.5395239384953d-5,2.5066282746310005d0/
      x=xx
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*dlog(tmp)-tmp
	ser=1.000000000190015d0
	DO j=1,6
		y=y+1.d0
		ser=ser+cof(j)/y
	END DO
	gammln=tmp+dlog(stp*ser/x)
      RETURN
      END      
            