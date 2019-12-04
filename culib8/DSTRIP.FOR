
c***********************************************************************
      SUBROUTINE DSTRIP(IQ,DR,K,F1,F2,FR,QNUM,OPT,KM,SL,C)
c
c     Calculates two nucleon transfer form factor via the
c     Bayman and Kallio method.
c***********************************************************************
c
c     IQ    Input quantum numbers for form factor
c     DR    Step size
c     K     Number of steps
c     F1    First  orbital
c     F2    Second orbital
c     FR    Output form factor
c     QNUM  Quantum numbers for orbitals
c     OPT
c
      implicit real*8(a-h,o-z)
      parameter (maxg=10, maxr=12)
      DIMENSION F1(400),F2(400),FR(800),QNUM(4,2),TVCC(10),IQ(3)
     1         ,D1(10),D2(10),C(32),AG(maxg),WG(maxg),BG(maxg)
     2         ,AR(maxr),WR(maxr)
C
      DATA KR,KX/0,0/
      data ag,wg,bg,ar,wr/maxg*0.,maxg*0.,maxg*0.,maxr*0.,maxr*0./
C
c     c1 = R1 scale
c     c2 = R2 scale
c     c3 = r1 scale
c     c4 = r2 scale
c     c5 = r0, integration scale
c     c6 = Pauli flag
c     c7 = order of integration
c     c8 = number of integration points
c
      IPAULI=C(6)
      R0    =C(5)
      L =IQ(1)
      IS=IQ(2)
      JJ=IQ(3)
      N1=QNUM(1,1)
      N2=QNUM(1,2)
      L1=QNUM(2,1)
      L2=QNUM(2,2)
      J1=QNUM(3,1)
      J2=QNUM(3,2)
      IS1=QNUM(4,1)
      IS2=QNUM(4,2)
      c3=c(3)
      c4=c(4)
      dr1=dr*c(1)
      dr2=dr*c(2)
      ITEMP=N1+N1+N2+N2
      NX=(((ITEMP+L1+L2+L+2)/4+2)+1+2)/2+2*IPAULI
      if(c(8).eq.0.0) then
        NR=(((ITEMP+L1+L2-L+3)/4+8)+1+2)/2+2*IPAULI
        NR=MIN0(NR,maxr)
        IF(NR.NE.KR) then
          CUT=0.0
          IERR=0
          ALFA=C(7)
          CALL GAUSSR(NR,KS,ALFA,AR,WR,IERR,CUT)
          NR=KS
          kr=nr
        endif
      else
        nr = c(8)
        do 45 i=1,nr
        ar(i)=c(i+8)
        wr(i)=c(i+nr+8)
        if(wr(i).ne.0.0) kr=i
   45   continue
      nr=kr
      endif
      NX=MIN0(NX,maxg)
      WRITE(6,9000)NR,NX
 9000 FORMAT(15H0 NO. R STEPS =,I3,18H     NO. X STEPS =,I3)
      IF(NX.NE.KX) then
        CALL LEGAUS(2*NX,AG,WG)
        KX=NX
        do 47 i=1,nx
        BG(i) = sqrt(1.0-ag(i)**2)
   47   continue
      endif
      LL=L+L
      LL1=L1+L1
      LL2=L2+L2
      FACT1=YXFCT(L1,LL1)/2.0**L1
      FACT2=YXFCT(L2,LL2)/2.0**L2
      FNORM=OPT*SQRT(FLOAT((LL1+1)*(LL2+1))/FLOAT(LL+1))
      IF(IPAULI.NE.0) GO TO 60
      TEMP=2.0
      DO 50 I=1,4
      IF(QNUM(I,1).NE.QNUM(I,2)) GO TO 55
   50 CONTINUE
      TEMP=4.0
   55 CONTINUE
      FNORM=FNORM/SQRT(TEMP)
   60 CONTINUE
      FNORM=FNORM*SQRT(FLOAT((LL+1)*(IS+1)*(J1+1)*(J2+1)))
     1     *WINEJ(LL1,IS1,J1,LL2,IS2,J2,LL,IS,JJ)
      FM1=1.0
      FM2=1.0
      FL1=L1
      FL2=L2
      LPL=MIN0(L1,L2)+1
      DO 80 M=1,LPL
      M2=M+M-2
      FM=M-1
      TVCC(M)=VCC(LL1,LL2,LL,M2,-M2)*2.0/SQRT(FM1*FM2)
      FM1=FM1*(FL1+FM+1)*(FL1-FM)
      FM2=FM2*(FL2+FM+1)*(FL2-FM)
   80 CONTINUE
      TVCC(1)=TVCC(1)/2.0
C
C     RS=r
C
C     R1=| C1*R+C3*r |
C     R2=| C2*R+C4*r |
C
C
C     CENTER OF MASS R LOOP
C
      R=0.0
      S=0.0
      DO 500 M=1,K
      R=R+DR1
      S=S+DR2
      RSQ=R**2
      SSQ=S**2
      SUMR=0.0
C
C     RELATIVE R LOOP
C
      DO 400 MR=1,KR
      RS=AR(MR)*R0
      SS=RS*C3
      RS=RS*C4
      RSSQ=RS**2+RSQ
      RPROD=2.0*R*RS
      SSSQ=SS**2+SSQ
      SPROD=2.0*S*SS
C
C     RELATIVE R ANGLE LOOP
C
      SUMX=0.0
      DO 300 MX=1,KX
      X =AG(MX)
      y =BG(mx)
      IX=0
  110 CONTINUE
      R1=SQRT(RSSQ+RPROD*X)
      R2=SQRT(SSSQ+SPROD*X)
      FK1=R1/DR
      K1=FK1
      K1=MAX0(K1,2)
      FK1=FK1-FLOAT(K1)
      IF(K1.GT.K) GO TO 300
      FK2=R2/DR
      K2=FK2
      K2=MAX0(K2,2)
      FK2=FK2-FLOAT(K2)
      IF(K2.GT.K) GO TO 300
      COS1=(R+RS*X)/R1
      COS2=(S+SS*X)/R2
      SIN1=abs(rs*y/r1)
      SIN2=abs(ss*y/r2)
  120 CONTINUE
      DM1=0.0
      DL1=FACT1
      IF(L1.EQ.0) GO TO 141
      DO 140 LI=1,L1
      DL1=DL1*SIN1
  140 CONTINUE
  141 CONTINUE
      D1(L1+1)=DL1
      DM2=0.0
      DL2=FACT2
      IF(L2.EQ.0) GO TO 151
      DO 150 LI=1,L2
      DL2=DL2*SIN2
  150 CONTINUE
  151 CONTINUE
      D2(L2+1)=DL2
      IF(L1.EQ.0) GO TO 171
      FJ1=1.0
      FL1=LL1
      FM1=LL1
      DO 170 LI=1,L1
      DK1=(FM1*COS1*DL1/SIN1-DM1)/(FJ1*FL1)
      FJ1=FJ1+1.0
      FL1=FL1-1.0
      FM1=FM1-2.0
      DM1=DL1
      DL1=DK1
      INDX=L1+1-LI
      D1(INDX)=DL1
  170 CONTINUE
  171 CONTINUE
      IF(L2.EQ.0) GO TO 181
      FJ2=1.0
      FL2=LL2
      FM2=LL2
      DO 180 LI=1,L2
      DK2=(FM2*COS2*DL2/SIN2-DM2)/(FJ2*FL2)
      FJ2=FJ2+1.0
      FL2=FL2-1.0
      FM2=FM2-2.0
      DM2=DL2
      DL2=DK2
      INDX=L2+1-LI
      D2(INDX)=DL2
  180 CONTINUE
  181 CONTINUE
      PROD=0.0
      DO 185 LI=1,LPL
      PROD=PROD+D1(LI)*D2(LI)*TVCC(LI)
  185 CONTINUE
  280 CONTINUE
C     INTERPOLATE - 4 POINT FORMULA
      FT1=-FK1*(FK1-1.)*(FK1-2.)*F1(K1-1)/6.
     1  +   (FK1**2-1.)*(FK1-2.)*F1(K1  )/2.
     2  -  FK1*(FK1+1.)*(FK1-2.)*F1(K1+1)/2.
     3  +        FK1*(FK1**2-1.)*F1(K1+2)/6.
      FT2=-FK2*(FK2-1.)*(FK2-2.)*F2(K2-1)/6.
     1  +   (FK2**2-1.)*(FK2-2.)*F2(K2  )/2.
     2  -  FK2*(FK2+1.)*(FK2-2.)*F2(K2+1)/2.
     3  +        FK2*(FK2**2-1.)*F2(K2+2)/6.
C     INTERPOLATE - 2 POINT FORMULA
C     FT1=FK1*(F1(K1+1)-F1(K1))+F1(K1)
C     FT2=FK2*(F2(K2+1)-F2(K2))+F2(K2)
c	if(m.eq.21) write(20,'(2i4,1p10e12.4)')mr,mx,x,r1,r2,r1**2+r2**2
c     1    ,ft1,ft2,ft1*ft2
      SUMX=SUMX+WG(MX)*PROD*FT1*FT2
      IF(IX.NE.0) GO TO 300
      IX=1
      IF(IPAULI.EQ.0) then
        ITEMP=K1
        K1=K2
        K2=ITEMP
        ATEMP=FK1
        FK1=FK2
        FK2=ATEMP
        ATEMP=COS1
        COS1=COS2
        COS2=ATEMP
        ATEMP=SIN1
        SIN1=SIN2
        SIN2=ATEMP
        IF(L1.EQ.L2) GO TO 280
        GO TO 120
      endif
      X=-X
      GO TO 110
  300 CONTINUE
      SUMR=SUMR+WR(MR)*SUMX
  400 CONTINUE
      SUMR=SUMR*FNORM
      FR(2*M-1)=FR(2*M-1)+SUMR
      IF(M.EQ.KM) SL=SUMR
  500 CONTINUE
 1000 CONTINUE
      RETURN
      END