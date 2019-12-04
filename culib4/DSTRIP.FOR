
      SUBROUTINE DSTRIP(IQ,DR,K,F1,F2,FR,QNUM,OPT,KM,SL,C)
c
c     Calculates two nucleon transfer form factor (Bayman and Kallio)
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
      parameter (maxg=10, maxr=12)
      DIMENSION F1(400),F2(400),FR(800),QNUM(4,2),IQ(4),TVCC(-9:9)
     1         ,D1(-9:9),D2(-9:9),C(32),AG(maxg),WG(maxg),BG(maxg)
     2         ,AR(maxr),WR(maxr)
C
      DATA KR,KX,krel/0,0,0/
      data ag, wg, bg, ar, wr/maxg*0.,maxg*0.,maxg*0.,maxr*0.,maxr*0./
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
      L    = IQ(1)
      IS   = IQ(2)
      JJ   = IQ(3)
      lrel = iq(4)
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
      NX=(((ITEMP+L1+L2+L+2)/4+2)+3+4*lrel)/2+2*IPAULI
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
      IF(NX .NE. KX .or. lrel .ne. krel) then
        CALL LEGAUS(2*NX,AG,WG)
        KX=NX
        krel=lrel
        DL1=sqrt(YXFCT(0,2*lrel))*YXFCT(lrel,2*lrel)/2.0**lrel
     1     *sqrt(float(2*lrel+1))*phasef(lrel)
        do 47 i=1,nx
        BG(i) = sqrt(1.0-ag(i)**2)
        DM1=0.0
        DO 46 LI=1,lrel
        DL1=DL1*BG(i)
   46 CONTINUE
        wg(i) = wg(i)*DL1
   47   continue
      endif
c      WRITE(20,'(15H0 NO. R STEPS =,I3,18H     NO. X STEPS =,I3)')NR,NX
c
      LL=L+L
      LL1=L1+L1
      LL2=L2+L2
      FACT1=YXFCT(L1,LL1)/2.0**L1
      FACT2=YXFCT(L2,LL2)/2.0**L2
      FNORM=OPT*SQRT(FLOAT((LL1+1)*(LL2+1))/FLOAT(LL+1))
     1     /vcc(LL,2*lrel,LL,0,2*lrel)
      IF(IPAULI. eq. 0) then
        TEMP=2.0
        DO 50 I=1,4
        IF(QNUM(I,1).NE.QNUM(I,2)) GO TO 55
   50   CONTINUE
        TEMP=4.0
   55   CONTINUE
        FNORM=FNORM/SQRT(TEMP)
      endif
c
      FNORM=FNORM*SQRT(FLOAT((LL+1)*(IS+1)*(J1+1)*(J2+1)))
     1     *WINEJ(LL1,IS1,J1,LL2,IS2,J2,LL,IS,JJ)
      if(l1 .le. l2) then
        mmin= -min0(l1,l2-lrel)
        mmax= l1
      else
        mmin=-min0(l2,l1-lrel)
        mmax= l2+lrel
      endif
      DO 80 M=mmin,mmax
      TVCC( M)=VCC(LL1,LL2,LL, 2*M,-2*M+2*lrel)
     1   *SQRT(yxfct(l1+abs( m     ),l1-abs( m     ))
     2        *yxfct(l2+abs(-m+lrel),l2-abs(-m+lrel)))
c	write(20,'(a,2i4,2f10.5)')' m1, m2, Tvcc(m1,m2)'
c     1	     ,m,-m+lrel,Tvcc( m)
   80 CONTINUE
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
      DO 500 M =1,K
      R=R+DR1
      S=S+DR2
      RSQ=R**2
      SSQ=S**2
      SUMR=0.0
C
C     RELATIVE R LOOP
C
c     RS = r/2
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
c
  120 CONTINUE
      DM1=0.0
      DL1=FACT1
      DO 140 LI=1,L1
      DL1=DL1*SIN1
  140 CONTINUE
c     time reversal phase
      D1( L1)=DL1*phasef(L1)
      D1(-L1)=DL1
      DM2=0.0
      DL2=FACT2
      DO 150 LI=1,L2
      DL2=DL2*SIN2
  150 CONTINUE
c     time reversal phase * exp(im*pi)
      D2( L2)=DL2
      D2(-L2)=DL2*phasef(L2)
      FJ1=1.0
      FL1=LL1
      FM1=LL1
      DO 170 m1=L1-1,0,-1
      DK1=(FM1*COS1*DL1/SIN1-DM1)/(FJ1*FL1)
      FJ1=FJ1+1.0
      FL1=FL1-1.0
      FM1=FM1-2.0
      DM1=DL1
      DL1=DK1
c     time reversal phase
      D1( m1)=DL1*phasef(m1)
      D1(-m1)=DL1
  170 CONTINUE
c	if(m .eq. 1 .and. mr .eq. 1) then
c	write(20,'(a,10(i4,f10.5))')' m1, D1 ='
c     1	     ,(m1,d1(m1),m1=l1,-l1,-1)
c	endif
      FJ2=1.0
      FL2=LL2
      FM2=LL2
      DO 180 m2=L2-1,0,-1
      DK2=(FM2*COS2*DL2/SIN2-DM2)/(FJ2*FL2)
      FJ2=FJ2+1.0
      FL2=FL2-1.0
      FM2=FM2-2.0
      DM2=DL2
      DL2=DK2
c     time reversal phase * exp(im*pi)
      D2( m2)=DL2
      D2(-m2)=DL2*phasef(m2)
  180 CONTINUE
c	if(m .eq. 1 .and. mr .eq. 1) then
c	write(20,'(a,10(i4,f10.5))')' m2, D2 ='
c     1	     ,(m2,d2(m2),m2=l2,-l2,-1)
c	endif
      PROD=0.0
      DO 185 LI=mmin,mmax
      PROD=PROD+D1( LI)*D2(-LI+lrel)*TVCC(LI)
c	if(m .eq. 1 .and. mr .eq. 1) then
c	write(20,'(a,2i4,4f10.5)')' m1, m2, D1, D2, Tvcc, Prod ='
c     1	     ,li,-li+lrel,d1(li),d2(-li+lrel),Tvcc(li),prod
c	endif
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
      SUMX=SUMX+WG(MX)*PROD*FT1*FT2
c	if(m.eq.34) write(20,'(2i4,1p10e12.4)')mr,mx,x,r1,r2,r1**2+r2**2
c     1    ,ft1,ft2,ft1*ft2,prod,wg(mx),sumx
      IF(IX.eq.0) then
        IX=1
        ITEMP=K1
        K1=K2
        K2=ITEMP
        ATEMP=FK1
        FK1=FK2
        FK2=ATEMP
        IF(L1.EQ.L2) GO TO 280
        ATEMP=COS1
        COS1=COS2
        COS2=ATEMP
        ATEMP=SIN1
        SIN1=SIN2
        SIN2=ATEMP
        GO TO 120
      endif
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