
      SUBROUTINE LEGAUS(LL,X,W)
C
C     GAUSS-LEGENDRE POINT AND WAIT ROUTINE FOR AN EVEN NO. OF POINTS
C
C	ll = order i.e. the total number of points for  -1< x <1
C	x  = points
C	w  = weights
C
      IMPLICIT REAL*8(A-H,O-Z)
      real*4 x,w
      DIMENSION X(100),W(100)
      Z3=-1.6/FLOAT(LL+1)
      Z2=3.0*Z3
      Z1=5.0*Z3
      NL=(LL+1)/2
      DO 200 L=1,NL
      ZOLD=0.0
      Z=Z1+3.0*(Z3-Z2)
      DO 50 J=1,10
      P1=0.0
      P2=1.0
      DO 30 I=1,LL
      P3=(FLOAT(I+I-1)*Z*P2-FLOAT(I-1)*P1)/FLOAT(I)
      P1=P2
      P2=P3
   30 CONTINUE
      DP=FLOAT(LL)*(P1-Z*P2)/(1.0-Z*Z)
      Z=Z-P2/DP
      IF(ABS(Z-ZOLD)/Z.LT.1.0E-10) GO TO 51
      ZOLD=Z
   50 CONTINUE
      WRITE(6,9100)L
 9100 FORMAT(28H0NO CONVERGENCE FOR ZERO NO.,I4)
   51 CONTINUE
      X(L)=Z
      W(L)=2.0/((1.0-Z*Z)*DP*DP)
      Z1=Z2
      Z2=Z3
      Z3=Z
  200 CONTINUE
      RETURN
      END