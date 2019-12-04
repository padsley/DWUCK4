
c***********************************************************************
      SUBROUTINE COU(R,RP,ETA,L,H,F,FP,G,GP,S)
c
c     Coulomb Function Subroutine
c
c     Calculates the functions at 2 points,  R and RP for convenience 
c     in matching to the output of some integration formulae.
c
c     R   First  argument
c     RP  Second argument
c     E   Coulomb parameter, eta
c     L   Number of angular momenta = lmax+1
c     H   Integration step size, usually = 0.10
c     F   Regular   function at R
c     FP  Regular   function at RP
c     G   Irregular function at R
c     GP  Irregular function at RP
c     S   Coulomb phase
c
c***********************************************************************
c
c      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(51),FP(51),G(51),GP(51),S(51)
C
      LL=L
      IF(LL.LT.50) THEN
        ELP=50.
        J=50
      ELSE
        ELP=LL
        J=LL
      ENDIF
      A=ATAN (ETA/ELP)
      B=SQRT (ETA**2+ELP**2)
      Y=A*(ELP-0.5)+ETA*( LOG(B)-1.)-SIN (A)/(12.*B)
     1  +SIN (3.*A)/(360.*B**3)-SIN (5.*A)/(1260.*B**5)
     2  +SIN (7.*A)/(1680.*B**7)-SIN (9.*A)/(1188.*B**9)
C
      K=J-1
      DO 100 I=1,K
      IF(J.LE.LL) S(J)=Y
      ELP=ELP-1.
      J=J-1
      Y=Y-ATAN (ETA/ELP)
  100 CONTINUE
      S(1)=Y
      DEL1=R-2.0*ETA
      RMAX=.41666667*(ETA**2+4.*ETA+3.)
      IF(RMAX.LT.10.0) RMAX=10.0
      DEL=R-RMAX
      IF(ETA.GE.5.0.AND.ABS (DEL1).LT.ABS (DEL)) THEN
      DEL=DEL1
        I1=2
C     Transition line expansion
      X=2.0*ETA
      T1=ETA**2
      T2=T1**2
      T3=ETA** .666666667
      T4=T3**2
      T5=T4**2
      T6=T3*T5
      T7=T4*T6
      T8=T3*T7
      T9=ETA** .166666667
      Y=1.22340402*T9*(1.+.495957017E-1/T4-.888888889E-2/T1
     1  +.245519918E-2/T6-.910895806E-3/T2+.845362E-3/T8)
      Z=-.707881773/T9*(1.-.172826039/T3+.317460317E-3/T1
     1  -.358121485E-2/T5+.311782468E-3/T2-.907396643E-3/T7)
      ELSE
C
      IF(DEL.GE.0.0.OR.ETA.EQ.0.0) THEN
        X=R
        I1=1
      ELSE
        X=RMAX
        I1=2
      ENDIF
C     Asymptotic expansion
      T1=ETA**2
      T2=2.*X
      T4=ETA/T2
      SS=1.
      TS=0.
      SL=0.
      TL=1.-ETA/X
      SSS=1.
      STS=0.
      SSL=0.
      STL=TL
      EN=0.
      DO 700 K=1,25
      IF(ABS (SS/SSS).GT.1.0E-10) THEN
        T5=EN+1.
        T6=T5+EN
        T7=EN*T5
        T8=T6*T4/T5
        T9=(T1-T7)/(T2*T5)
        T5=T8*SS-T9*TS
        TS=T8*TS+T9*SS
        SS=T5
        T5=T8*SL-T9*TL-SS/X
        TL=T8*TL+T9*SL-TS/X
        SL=T5
        SSS=SSS+SS
        STS=STS+TS
        SSL=SSL+SL
        STL=STL+TL
        EN=EN+1.
      ENDIF
  700 CONTINUE
C
      T3=X-ETA* LOG(T2)+S(1)
      T8=SIN (T3)
      T9=COS (T3)
      Y=SSS*T9-STS*T8
      Z=SSL*T9-STL*T8
      ENDIF
C
      DO 810 I=1,I1
      IF(I.EQ.I1) THEN
        G(1) = Y
        W    = Z
        DEL = RP - R
      ENDIF
C     Runge-Kutta integration
      N=ABS (DEL/H)+1.0
      DX=DEL/FLOAT(N)
      T1=DX/2.
      T2=DX/8.
      T3=2.0*ETA
      DO 805 J=1,N
      T4=DX*(T3/X-1.)*Y
      X=X+T1
      T5=DX*(T3/X-1.)*(Y+T1*Z+T2*T4)
      X=X+T1
      T6=DX*(T3/X-1.)*(Y+DX*Z+T1*T5)
      Y=Y+DX*(Z+(T4+2.*T5)/6.)
      Z=Z+(T4+4.*T5+T6)/6.
  805 CONTINUE
  810 CONTINUE
      GP(1)=Y
C
      T1=ETA**2
      T8=SQRT (1.+T1)
      G(2)=((1./R+ETA)*G(1)-W)/T8
      GP(2)=((1./RP+ETA)*Y-Z)/T8
      T2=1.
      T3=2.
C     Recur irregular function upwards
      DO 910 I=3,LL
      T4=T2+T3
      T5=T2*T3
      T6=T3*SQRT (T2**2+T1)
      T7=T2*SQRT (T3**2+T1)
      G (I)=(T4*(ETA+T5/R )*G (I-1)-T6*G (I-2))/T7
      GP(I)=(T4*(ETA+T5/RP)*GP(I-1)-T6*GP(I-2))/T7
      T2=T2+1.
      T3=T3+1.
  910 CONTINUE
      N=MAX0(L+11, INT(2.0*R+11.0) )
      I1=N
      Y =1.0E-20
      YP=Y
      Z =0.
      ZP=Z
C     Recur regular function downwards
      DO 930 I=1,I1
      T2=N
      T3=T2+1.
      T4=T2+T3
      T5=T2*T3
      T6=T2*SQRT (T3**2+T1)
      T7=T3*SQRT (T2**2+T1)
      X =(T4*(ETA+T5/R )*Y -T6*Z )/T7
      XP=(T4*(ETA+T5/RP)*YP-T6*ZP)/T7
      IF(N.LE.LL) THEN
        F (N)=X
        FP(N)=XP
      ELSE
c     Scale for l > lmax
        IF(ABS (X).GT.1.0) THEN
          Y =Y *1.0E-20
          YP=YP*1.0E-20
          X =X *1.0E-20
          XP=XP*1.0E-20
        ENDIF
      ENDIF
      N=N-1
      Z =Y
      ZP=YP
      Y =X
      YP=XP
  930 CONTINUE
      Z =1./(T8*(F (1)*G (2)-F (2)*G (1)))
      ZP=1./(T8*(FP(1)*GP(2)-FP(2)*GP(1)))
      DO 950 I=1,LL
      F (I)=F (I)*Z
  950 FP(I)=FP(I)*ZP
      RETURN
       END