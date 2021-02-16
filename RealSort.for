      SUBROUTINE RealSort(n,a)
C
C From Numerical Recipes in Fortran (called "shell")
C
      IMPLICIT NONE
C
      INTEGER*4 n
      INTEGER*4 i,j,inc
      REAL*8 v,a(n)
C
C
      inc=1
1     inc=3*inc+1
      if(inc.le.n)goto 1
2     continue
        inc=inc/3
        do 11 i=inc+1,n
          v=a(i)
          j=i
3         if(a(j-inc).gt.v)then
            a(j)=a(j-inc)
            j=j-inc
            if(j.le.inc)goto 4
          goto 3
          endif
4         a(j)=v
11      continue
      if(inc.gt.1)goto 2
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $2#!'U,)'=.
