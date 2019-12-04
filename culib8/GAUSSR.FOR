
c***********************************************************************
      SUBROUTINE GAUSSR(NMAX,INDEX,ALFA,AG,WG,IERR,CUTOFF)
C
c     Gauss-Hermite and gauss-laguerre point and weight routine
c***********************************************************************
C
C     IF ALFA IS INTEGER       -- GAUSS-LAGUERRE
C     IF ALFA IS INTEGER + 1/2 -- GAUSS-HERMITE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AG(100),WG(100)
      data eps/1.e-6/
      DATA PI,SQRPI/3.14159265,1.77245385/
c
      INDEX=0
      FI=NMAX
      FKI=4.0*(FI+(ALFA+1.0)*0.5)
      FLN= LOG(FI)
C
      K=ALFA
      JJ=(ALFA-FLOAT(K))*2.0
      DY=0.0
      IF(JJ.NE.0) GO TO 19
      FNORM=1.0
      IF(K.EQ.0) GO TO 11
      DO 10 J=1,K
      FNORM=FNORM*(ALFA+FI+1.0-FLOAT(J))
   10 CONTINUE
   11 CONTINUE
      GO TO 25
   19 CONTINUE
      FNORM=SQRPI/2.0
      DO 20 J=1,NMAX
      FNORM=FNORM*(1.0+ALFA/FLOAT(J))
   20 CONTINUE
      K=ALFA+1.0
      IF(K.EQ.0) GO TO 25
      DO 22 J=1,K
      FNORM=FNORM*(FLOAT(J)-0.5)
   22 CONTINUE
   25 CONTINUE
      Y=0.0
      Z1=0.0
      Z2=0.0
C
      DO 80 J=1,NMAX
      FJ=J
      FKJ=4.0*(FJ+(ALFA+1.0)*0.5)
      Z=((FJ+0.5*ALFA-0.25)*PI)**2/FKI
      Z=Z*(1.0+Z*(1.0+PI* LOG(FJ)*FLN/(8.0*(FI+FLN-FJ+eps)))/(3.0*FKI))       
      DELZ=Z-Z1
      Z1=Z
      Z=Y+DELZ
      Z3=Z
      DO 74 M=1,20
      A1=0.0
      A2=1.0
      DO 70 K=1,NMAX
      FK=K
      A3=((2.0*FK+ALFA-1.0-Z)*A2-(FK+ALFA-1.0)*A1)/FK
      A1=A2
      A2=A3
      B2=(FK*A2-(FK+ALFA)*A1)/Z
   70 CONTINUE
      Y=Z-A2/B2
      IF(ABS((Z-Y)/Z).LT.3.E-14) GO TO 75
      Z=Y
   74 CONTINUE
      IERR=J
   75 CONTINUE
      DZ=Z-Z2
      Z2=Z
      DZ=((FNORM/Z)/B2)/B2
      IF(DZ.LT.CUTOFF.AND.DZ.LT.DY) GO TO 100
      DY=DZ
      INDEX=J
      AG(J)=Y
      WG(J)=DZ
      IF(JJ.NE.0) AG(J)=SQRT(Y)
   80 CONTINUE
  100 CONTINUE
      RETURN
      END