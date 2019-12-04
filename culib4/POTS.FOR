c***********************************************************************
      SUBROUTINE POTS(U,V)
c
c     Calculates the potentials or form factors
c***********************************************************************
c
c     IMPLICIT REAL*8(A-H,O-Z)
      COMMON/POTTER/DRX,AFACT(2),VFACT,SFACT,ENG,RM,G(4),ETAX,FKX,ETAKX
     1             ,RCX,HBARC2,ABETA(3),FLDF(3)
     2             ,NX,LAMX,KMXX,KX,IBX,LPLUSX,ICON4,NSPC,IDIRAC,ICHK
 
      CHARACTER*18 B(26)
      DIMENSION UT(5),CN(16),CP(16),YLAM(16),PLM(8)
     1         ,XG(8),WG(8),U(800),V(800),LDFRM(3)
C      DIMENSION UD(800)
      EQUIVALENCE (YLAM(1),PLM(1))
C
c
      DATA b /                ' NX=0 No potential', ' NX=1 VOLUME W-S  '
     1, ' NX=2 SURFAC, W-S ', ' NX=3 2ND DERIV   ', ' NX=4 L.S VOLUME  '
     2, ' NX=5 L.S SURFACE ', ' NX=6 VOL*R**POWR ', ' NX=7 SUR*R**POWR '
     3, ' NX=8 EXTERN FORMF', ' NX=9 HARMONIC OSC', ' NX=10 GAUSSIAN   '
     4, ' NX=11 DEFORM VOL ', ' NX=12 DEFORM SURF', ' NX=13 Hulthen    '
     5, ' NX=14 Yukawa Lam ', ' NX=15 Yukawa L*S ', ' NX=16 NO OPTION  '
     6, ' NX=17 NO OPTION  ', ' NX=18 NO OPTION  ', ' NX=19 NO OPTION  '
     7, ' NX=20 VECTOR     ', ' NX=21 SCALAR     ', ' NX=22 Not used   '
     8, ' NX=23 Not used   ', ' NX=24 Not used   ', ' NX=25 Not used   '
     9/
C
      DATA NG,NGX/8,0/
C
      ETA4 = 6.0
      ETA5 =10.0
      SQRPI= 1.772453851
      PI=3.141592654
C
      IVFLAG=0
      ISFLAG=0
      FACT=VFACT
   70 CONTINUE
C
C     READ IN CARD SET 5,6,OR 7   POTENTIAL CARDS
C
      READ (5,9000)FZ,VR,RY,AR,VSOR,VI,RZ,AI,VSOI,PWR
C
      NZ=ABS(FZ)
      RR=ABS(RY)*AFACT(1)
      IF(RY.LT.0.0) RR=RR+ABS(RY)*AFACT(2)
      RI=ABS(RZ)*AFACT(1)
      IF(RZ.LT.0.0) RI=RI+ABS(RZ)*AFACT(2)
      IF(ICON4.NE.2.OR.NSPC.LT.3) THEN
        WRITE(6,9509)B(NZ+1),VR,RY,AR,RR,VSOR
        WRITE(6,9510)        VI,RZ,AI,RI,VSOI,PWR
      ENDIF
C
      KFLAG=0
      IF(NX.LT.3) THEN
        VR=VR*FACT
        VI=VI*FACT
        KT=FKX  * MAX (RR,RI)+ETA5
        LPLUSX=MAX0(LPLUSX,KT)
        KT=(2.3*ETA4* MAX (AR,AI)+ MAX (RR,RI))/DRX
      ELSE
        IF(ENG.EQ.0.0) THEN
          KT=(2.3*ETA4* MAX (AR,AI)+ MAX (RR,RI))/DRX
        ELSE
          RM= MAX (RM,RR)
          RM= MAX (RM,RI)
          IF(RM.EQ.0.0) RM=1.0
          VR=VR*FACT
          VI=VI*FACT
          KT=(2.3*ETA4/SQRT(FKX**2+2.0*ETAX*FKX/RM))/DRX
        ENDIF
      ENDIF
      KX=MIN0(MAX0(KX,KT),KMXX)
   83 CONTINUE
      IF(AR.EQ.0.0) THEN
        F1=0.0
        F2=0.0
      ELSE
        F2=EXP(-DRX  /AR)
        F1=EXP( RR/AR)
      ENDIF
      IF(AI.EQ.0.0) THEN
        F3=0.0
        F4=0.0
      ELSE
        F4=EXP(-DRX  /AI)
        F3=EXP( RI/AI)
      ENDIF
C
      IF(NX.GE.3.AND.ENG.EQ.0.0.AND.NZ.LE.5) THEN
        IF(AR.NE.0.0) VR=VR*(RR/AR)**(PWR+1)
        IF(AI.NE.0.0) VI=VI*(RI/AI)**(PWR+1)
      ENDIF
      IF(NZ.EQ.0) GO TO 6050
      GO TO ( 100, 200, 300, 400, 500, 600, 700, 800, 900,1000),NZ
      GO TO (1100,1200,1300,1400,1500,1600,1700,1800,1900,2000),NZ-10
      GO TO (2100),NZ-20
      write(6,'(''0Invalid potential option '',i3)')nz
      
C
C     VOLUME WOODS SAXON
C
  100 CONTINUE
      DO 160 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3/(1.0+F3)
  160 CONTINUE
      GO TO 6000
C
C     1ST DERIVATIVE WOODS SAXON
C
  200 CONTINUE
      DO 260 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )+VR*F1/(1.0+F1)**2
      U(MK+1)=U(MK+1)+VI*F3/(1.0+F3)**2
  260 CONTINUE
      GO TO 6000
C
C     2ND DERIVATIVE WOODS SAXON
C
  300 CONTINUE
      DO 360 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1*(1.0-F1)/(1.0+F1)**3
      U(MK+1)=U(MK+1)-VI*F3*(1.0-F3)/(1.0+F3)**3
  360 CONTINUE
      GO TO 6000
C
C     L.S VOLUME WOODS SAXON
C
  400 CONTINUE
      IBX=1
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      R=0.0
      DO 460 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1/(R*(1.0+F1)**2)
      V(MK+1)=V(MK+1)-VI*F3/(R*(1.0+F3)**2)
  460 CONTINUE
      GO TO 6000
C
C     L.S 1ST DERIVATIVE WOODS SAXON
C
  500 CONTINUE
      IBX=1
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      R=0.0
      DO 560 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )+VR*F1*(1.0-F1)/(R*(1.0+F1)**3)
      V(MK+1)=V(MK+1)+VI*F3*(1.0-F3)/(R*(1.0+F3)**3)
  560 CONTINUE
      GO TO 6000
C
C     WOOD SAXON*R**PWR
C
  600 CONTINUE
      R=0.0
      DO 660 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      RPWR=R**PWR
      U(MK  )=U(MK  )-VR*F1*RPWR/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3*RPWR/(1.0+F3)
  660 CONTINUE
      GO TO 6000
C
C      1ST DERIVATIVE WOOD SAXON*R**PWR
C
  700 CONTINUE
      R=0.0
      DO 760 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      RPWR=R**PWR
      U(MK  )=U(MK  )+VR*F1*RPWR/(1.0+F1)**2
      U(MK+1)=U(MK+1)+VI*F3*RPWR/(1.0+F3)**2
  760 CONTINUE
      GO TO 6000
C
C     EXTERNAL FORM FACTOR
C
  800 CONTINUE
      IF(NX.GE.3) THEN
        READ (5,9000)G
        WRITE(6,9508)G
      ENDIF
      READ(5,9000)F1,F2
C
C     F2 = 0  REAL CENTRAL
C     F2 = 1  IMAG CENTRAL
C     F2 = 2  REAL SPIN ORBIT
C     F2 = 3  IMAG SPIN ORBIT
C
      IF(F2.EQ.0.0.OR.F2.EQ.2) THEN
        F3=VR
        MK=1
      ELSE
        F3=VI
        MK=2
      ENDIF
      IF(F3.EQ.0.0) F3=1.0
      KT=F1
      DO 820 M=1,KT,5
      READ(5,9100)UT
      IF(M.GT.KMXX) GO TO 820
      DO 810 I=1,5
      IF(F2.LT.2.0) THEN
        U(MK  )=U(MK  )+UT(I)*F3
      ELSE
        V(MK  )=V(MK  )+UT(I)*F3
      ENDIF
      MK=MK+2
  810 CONTINUE
  820 CONTINUE
C
      KX=min0(KT,kmxx)
 
      GO TO 6000
C
C     HARMONIC OSCILLATOR, NORMALIZED   N*EXP(-(R/RY)**2/2)
C
  900 CONTINUE
      READ (5,9000)G
      WRITE(6,9508)G
      F1=1.0/RY**2
      F2=F1/RY
      F3=0.5
      F4=SQRPI*0.5
      L=G(2)
      IF(L.NE.0) THEN
        DO 930 I=1,L
        F3=F3+1.0
        F4=F4*F3
        F2=F2*F1
  930   CONTINUE
      ENDIF
      NN=G(1)
      T1=1.0
      T2=F4
C     LAGUERRE POLYNOMIAL COEFFICIENTS = (Abramowitz-Stegun)*(-1)**nn
      CN(1)=(-1.0)**NN
      IF(NN.NE.0) THEN
        DO 940 I=1,NN
        F3=F3+1.0
        T1=T1*FLOAT(I)
        T2=T2*F3
        CN(I+1)=-CN(I)*F1*FLOAT(NN+1-I)/(FLOAT(I)*F3)
  940   CONTINUE
      ENDIF
      ANORM=SQRT(2.0*F2*T1/T2)*T2/(T1*F4)
      IF(VR.NE.0.0) ANORM=VR*ANORM
      KT=10.0*RY/DRX
      KT=MIN0(KT,KMXX)
      R=0.0
      F1=F1/2.0
      DO 960 M=1,KT
      MK=M+M-1
      R=R+DRX
      R2=R**2
      F2=CN(1)
      IF(NN.EQ.0) GO TO 951
      F3=1.0
      DO 950 I=1,NN
      F3=F3*R2
      F2=F2+CN(I+1)*F3
  950 CONTINUE
  951 CONTINUE
      U(MK  )=U(MK  )+F2*ANORM*EXP(-F1*R2)*R**L
  960 CONTINUE
      GO TO 6000
C
C     EXP(-(R/R0)**2)*R**POWR
C
 1000 CONTINUE
      IF(VR.NE.0.0) THEN
        R=0.0
        DO 1060 M=1,KX
        MK=M+M-1
        R=R+DRX
        U(MK  )=U(MK  )-VR*EXP(-(R/RY)**2)*R**PWR
 1060 CONTINUE
      endif
      IF(VI.NE.0.0) THEN
        R=0.0
        DO 1070 M=1,KX
        MK=M+M-1
        R=R+DRX
        U(MK+1)=U(MK+1)-VI*EXP(-(R/RZ)**2)*R**PWR
 1070 CONTINUE
      endif
      GO TO 6000
C
C     DEFORMED VOLUME OR SURFACE BY YLM EXPANSION
C
 1100 CONTINUE
      IF(NGX.NE.NG) THEN
        CALL LEGAUS(2*NG,XG,WG)
        NGX=NG
      ENDIF
      T2=(-1.0)**LAMX
      IF(ICHK.EQ.0) THEN
        READ (5,9000) (ABETA(J),FLDF(J) ,J=1,3)
        WRITE(6,9000)
        WRITE(6,9512) (ABETA(J),FLDF(J) ,J=1,3)
      ENDIF
      LMAX=LAMX+1
      DO 1101 J=1,3
      LDFRM(J)=FLDF(J)
      LMAX=MAX0(LMAX,LDFRM(J)+1)
 1101 CONTINUE
      T2=(-1.0)**LAMX
C
      DO 1108 I=1,NG
      CN(I  )=0.0
      CN(I+NG)=0.0
      P2=0.0
      P3=1.0
C
      DO 1106 M=1,LMAX
      L=M-1
      FL=L-1
      IF(L.EQ.0) GO TO 1102
      P3=((2.0*FL+1.0)*XG(I)*P2-FL*P1)/(FL+1.0)
 1102 CONTINUE
      DO 1103 J=1,3
      IF(ABETA(J).EQ.0.0) GO TO 1103
      IF(L.NE.LDFRM(J)) GO TO 1103
      FACTOR=P3*ABETA(J)*SQRT(FLOAT(L+L+1)/(4.0*PI))
      CN(I  )=CN(I  )+FACTOR
      CN(I+NG)=CN(I+8)+FACTOR*(-1.0)**LDFRM(J)
 1103 CONTINUE
      IF(L.NE.LAMX) GO TO 1104
      YLAM(I  )= P3*WG(I)*SQRT(FLOAT(L+L+1)*PI)
      YLAM(I+NG)=YLAM(I  )*T2
 1104 CONTINUE
      P1=P2
      P2=P3
 1106 CONTINUE
 1108 CONTINUE
      DO 1170 I=1,2
      IF(I.EQ.1) THEN
        IF(VR.EQ.0.0) GO TO 1170
        VX=VR
        RX=RR
        AX=AR
        F1=1.0
        IFL=-1
      ELSE
        IF(VI.EQ.0.0) GO TO 1170
        VX=VI
        RX=RI
        AX=AI
        F1=1.0
        F2=F4
        IFL=0
      ENDIF
      DO 1135 J=1,16
      CP(J)=EXP((1.0+CN(J))*RX/AX)
 1135 CONTINUE
      IF(LAMX.EQ.0) VX=VX/(SQRPI*2.0)
      J=NZ-10
      IF(J.EQ.2) GO TO 1151
      DO 1150 M=1,KX
      MK=M+M+IFL
      VTEMP=0.0
      F1=F1*F2
      DO 1145 J=1,16
      F3=F1*CP(J)
      VTEMP=VTEMP-YLAM(J)*VX*F3/(1.0+F3)
 1145 CONTINUE
      U(MK  )=U(MK  )+VTEMP
 1150 CONTINUE
      GO TO 1170
 1151 CONTINUE
      DO 1160 M=1,KX
      MK=M+M+IFL
      VTEMP=0.0
      F1=F1*F2
      DO 1155 J=1,16
      F3=F1*CP(J)
      VTEMP=VTEMP+YLAM(J)*VX*F3/(1.0+F3)**2
 1155 CONTINUE
      U(MK  )=U(MK  )+VTEMP
 1160 CONTINUE
 1170 CONTINUE
      GO TO 6000
 1200 CONTINUE
      GO TO 1100
C
C     HULTHEN WAVE FUNCTION= NORM*(EXP(-R/RY)-EXP(-R/RZ))/R
C
 1300 CONTINUE
      READ (5,9000)G
      WRITE(6,9508)G
      T1=1.0/RY
      T2=1.0/RZ
      T3=T2**2-T1**2
      T4=SQRT(2.0*T1*(T1+T2)*T2)/(T2-T1)
      IF(VR.NE.0.0) KT=16.0* MIN (RY,RZ)/DRX
      IF(VR.EQ.0.0) KT=16.0* MAX (RY,RZ)/DRX
      KX=MIN0(KX,KMXX)
      F1=1.0
      F2=EXP(-DRX  *T1)
      F3=1.0
      F4=EXP(-DRX  *T2)
      R=0.0
      DO 1260 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      V(M)=T3*F3/(F1-F3)
      TEMP=1.0
      IF(VR.NE.0.0) TEMP=V(M)/FACT
      U(MK  )=U(MK  )+TEMP*T4*(F1-F3)/R
 1260 CONTINUE
      GO TO 6000
c
c     Yukawa L=LAM
c
c     v = V-0*exp(-r/rx)/(r/rx)   ay < r 
c       = Wood-Saxon              r  < ay
c
 1400 CONTINUE
	lam=lamx
      if(ry.eq.0.0) ry=1.0
      if(rz.eq.0.0) rz=1.0
      if(ar.eq.0.0) ar=drx
      if(ai.eq.0.0) ai=drx
      f1=1.0
      f2=exp(-drx/ry)
      xr=ar/ry
      t1=1.0+1.0/xr
      t2=1.0
c     Recur Hankel functions
      do 1420 i=0,lam
      t0=t1
      t1=t2
      t2=float(2*i-1)*t1/xr+t0
 1420 continue
      vzr  =2.0*vr*t2*exp(-xr)/xr
      betar=2.0*(float(lam+1)*t2/xr+t1)/t2
      f3=1.0
      f4=exp(-drx/rz)
      xi=ai/rz
      t1=1.0+1.0/xi
      t2=1.0
c     Recur Hankel functions
      do 1430 i=0,lam
      t0=t1
      t1=t2
      t2=float(2*i-1)*t1/xi+t0
 1430 continue
      vzi  =2.0*vi*t2*exp(-xi)/xi
      betai=2.0*(float(lam+1)*t2/xi+t1)/t2
c
      if(vr.ne.0.0) then
      r =0.0
      do 1460 m=1,k
      r=r+drx
      mk=m+m-1
      f1=f1*f2
      f3=f3*f4
      if(r.lt.ar) then
        u(mk  )=u(mk  )-vzr/(1.0+exp(betar*(r/ry-xr)))
      else
      xr=r/ry
      t1=1.0+1.0/xr
      t2=1.0
c     Recur Hankel functions
      do 1450 i=0,lam
      t0=t1
      t1=t2
      t2=float(2*i-1)*t1/xr+t0
 1450 continue
        u(mk  )=u(mk  )-vr*t2*f1/xr
      endif
 1460 continue
      endif
      if(vi.ne.0.0) then
      r =0.0
      do 1480 m=1,k
      r=r+drx
      mk=m+m-1
      if(r.lt.ai) then
        u(mk+1)=u(mk+1)-vzi/(1.0+exp(betai*(r/rz-xi)))
      else
      xi=r/rz
      t1=1.0+1.0/xi
      t2=1.0
c     Recur Hankel functions
      do 1470 i=0,lam
      t0=t1
      t1=t2
      t2=float(2*i-1)*t1/xi+t0
 1470 continue
        u(mk+1)=u(mk+1)-vi*t2*f3/xi
      endif
 1480 continue
      endif
      go to 6000
c
c     Yukawa L*S
c
 1500 CONTINUE
      r =0.0
      if(ry.eq.0.0) ry=1.0
      if(rz.eq.0.0) rz=1.0
      if(ar.eq.0.0) ar=1.0
      if(ai.eq.0.0) ai=1.0
      f1=1.0
      f2=exp(-drx/ry)
      xr=ar/ry
      betar=2.0*(1.0+3.0/xr+3.0/xr**2)/(1.0+1.0/xr)
      vzr  =2.0*vr*exp(-xr)*(1.0+1.0/xr)/xr**2
      f3=1.0
      f4=exp(-drx/rz)
      xi=ai/rz
      betai=2.0*(1.0+3.0/xi+3.0/xi**2)/(1.0+1.0/xi)
      vzi  =2.0*vi*exp(-xi)*(1.0+1.0/xi)/xi**2
      do 1560 m=1,k
      r=r+drx
      mk=m+m-1
      f1=f1*f2
      f3=f3*f4
      if(r.lt.ar) then
        u(mk  )=u(mk  )-vzr/(1.0+exp(betar*(r/ry-xr)))
      else
        u(mk  )=u(mk  )-vr *f1*(1.0+1.0*ry/r)/(r/ry)**2
      endif
      if(r.lt.ai) then
        u(mk+1)=u(mk+1)-vzi/(1.0+exp(betai*(r/rz-xi)))
      else
        u(mk+1)=u(mk+1)-vi *f3*(1.0+1.0*rz/r)/(r/rz)**2
      endif
 1560 continue
      go to 6000
 1600 CONTINUE
 1700 CONTINUE
       go to 6000
 1800 CONTINUE
      DO 1860 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1/(1.0+F1)
      V(MK+1)=V(MK+1)-VI*F3/(1.0+F3)
 1860 CONTINUE
      GO TO 6000
 1900 CONTINUE
c
c     Coulomb excitation for a deformed uniform charge distribution
c
      IF(ICHK.EQ.0) THEN
        READ (5,9000) (ABETA(J),FLDF(J) ,J=1,3)
        WRITE(6,9000)
        WRITE(6,9512) (ABETA(J),FLDF(J) ,J=1,3)
      ENDIF
      if(vr.eq.0.0) vr=1.0
      do 1990 i=1,3
      beta=abeta(i)
      if(beta.ne.0.0) then
c     set flag for unbound stripping evaluation of coulomb ex.
      ibx=4
      flam1=beta*sqrt(float(2*lamx+1)/(4.0*pi))
      flam2=flam1*phasef(lamx)
c
      if(ngx.ne.ng) then
        ngx=ng
        call legaus(2*ng,xg,wg)
      endif
c
      an=0.0
      do 1920 k=1,ng
      p1=0.0
      p2=1.0
      if(lamx.gt.0) then
        do 1910 j=1,lamx
        p3=(float(2*j-1)*xg(k)*p2 - float(j-1)*p1)/float(j)
        p1=p2
        p2=p3
 1910 continue
      endif
      plm(k)=p2
c
c     calculate normalization
c
      r1=rcx*(1.0+flam1*p2)
      r2=rcx*(1.0+flam2*p2)
      an=an+(r1**3+r2**3)*wg(k)
 1920 continue
c
      an =an/3.0
c
      r  =0.0
c
      do 1940 m=1,kx
      mk=m+m-1
      r=r+drx
      sum=0.0
      do 1930 k=1,ng
      r1=rcx*(1.0+flam1*plm(k))
      r2=rcx*(1.0+flam2*plm(k))
c
      if(r.gt.r1) then
        s1=r1**(lamx+3)/(r**(lamx+1)*float(lamx+3))
      else
        if(lamx.eq.2) then
          s1=( log(r1/r)+1.0/float(lamx+3))*r**2
        else
          s1=(r**lamx/r1**(lamx-2)
     1      -float(2*lamx+1)*r**2/float(lamx+3))/float(2-lamx)
        endif
      endif
c
      if(r.gt.r2) then
        s2=r2**(lamx+3)/(r**(lamx+1)*float(lamx+3))
      else
        if(lamx.eq.2) then
          s2=( log(r2/r)+1.0/float(lamx+3))*r**2
        else
          s2=(r**lamx/r2**(lamx-2)
     1      -float(2*lamx+1)*r**2/float(lamx+3))/float(2-lamx)
        endif
      endif
c
      sum=sum+(s1+s2*phasef(lamx))*wg(k)*plm(k)
 1930 continue
c
      sum=sum/an
      if(beta.ne.0.0.and.lamx.ne.0) sum=sum/flam1
      u(mk  )=u(mk  )+sum*(vr*etakx/fact)
 1940 continue
      endif
 1990 continue
      GO TO 6000
C
C     VECTOR POTENTIAL
C     VSOR, VSOI ARE THE THIRD PARAMETERS
C     IN THE 3 PARAMETER FERMI MODEL
C         [1 + VSO?*(R/R?)**2]
C
 2000 CONTINUE
      IVFLAG=-1
      R=0.0
      DO 2060 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1*(1.0+VSOR*(R/RR)**2)/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3*(1.0+VSOI*(R/RI)**2)/(1.0+F3)
 2060 CONTINUE
      GO TO 6020
C
C     SCALAR POTENTIAL
C     VSOR, VSOI ARE THE FERMI THIRD PARAMETERS
C     IN THE 3 PARAMETER FERMI MODEL
C         [1 + VSO?*(R/R?)**2]
C
 2100 CONTINUE
      ISFLAG=-1
      VR=VR*(SFACT/VFACT)
      VI=VI*(SFACT/VFACT)
      R=0.0
      DO 2160 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1*(1.0+VSOR*(R/RR)**2)/(1.0+F1)
      V(MK+1)=V(MK+1)-VI*F3*(1.0+VSOI*(R/RI)**2)/(1.0+F3)
 2160 CONTINUE
      GO TO 6020
C
C      END OF POTENTIALS
C
 6000 CONTINUE
      IDIRAC=1
 6020 CONTINUE
      IF(KFLAG.NE.0.OR.NZ.GT.5) GO TO 6050
      IF(ABS(VSOR)+ABS(VSOI).EQ.0.0) GO TO 6050
      NZ=NZ+3
      KFLAG=1
      VR=VR*VSOR/45.2
      VI=VI*VSOI/45.2
      GO TO 83
 6050 CONTINUE
      IF(FZ.GT.0.0) GO TO 70
C
C     PROCESS DIRAC POTENTIALS
C
C     ENTRY WITH  U -> K**2 - VFACT*V   VFACT = 2.0*W1 /HBARC**2
C                 V ->      - SFACT*S   SFACT = 2.0*FM1/HBARC**2
C
      IF(IVFLAG.NE.0.AND.ISFLAG.NE.0) THEN
      IF(IDIRAC.EQ.1) THEN
        WRITE(6,9515)
      ENDIF
      IDIRAC=-1
C
C      KT2=KX+KX
C      WRITE(20,7777)' ENTRY  POTENTIALS'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)
C      WRITE(20,7778)(V(I),I=1,KT2)
      W1M1=(VFACT+SFACT)*0.5*HBARC2
 
      DO 6100 M=1,Kx
      MK=M+M-1
      VVR=(U(MK  )-FKX**2)/VFACT
      VVI= U(MK+1)        /VFACT
      VSR= V(MK  )        /SFACT
      VSI= V(MK+1)        /SFACT
      U(MK  )=U(MK  )+V(MK  ) +(VVR**2-VVI**2 - VSR**2+VSI**2)/HBARC2
      U(MK+1)=U(MK+1)+V(MK+1) +(VVR*VVI - VSR*VSI)*2.0/HBARC2
      T1 = W1M1 + VVR-VSR
      T2 =        VVI-VSI
      V(MK  )=0.5* LOG (T1**2 + T2**2)
      V(MK+1)=ATAN2(T2,T1)
 6100 CONTINUE
C      WRITE(20,7777)' SECOND POTENTIALS'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)
C      WRITE(20,7778)(V(I),I=1,KT2)
C
      R=FLOAT(KX+1)*DRX
      MK=KX+KX-1
      D2=V(MK  )
      D1=D2
      A2=V(MK+1)
      A1=A2
      IBX=1
      DO 6150 M=2,KX
      MK=MK-2
      R=R-DRX
      D3=D2
      D2=D1
      D1=V(MK  )
      A3=A2
      A2=A1
      A1=V(MK+1)
C     FIRST  DERIVATIVE TERMS
      DPR=(D3-D1)/(2.0*DRX)
      DPI=(A3-A1)/(2.0*DRX)
      V(MK+2)=2.0*DPR/R
      V(MK+3)=2.0*DPI/R
C     SECOND DERIVATIVE TERMS
      DPPR=(D3+D1-2.0*D2)/DRX**2
      DPPI=(A3+A1-2.0*A2)/DRX**2
      UDR =0.5*DPPR-0.25*(DPR**2-DPI**2)+DPR/R
      UDI =0.5*DPPI-0.25*(2.0*DPR*DPI  )+DPI/R
C      UD(MK+2)=UDR
C      UD(MK+3)=UDI
      U(MK+2)=U(MK+2)+UDR
      U(MK+3)=U(MK+3)+UDI
 6150 CONTINUE
      V(1   )=V(3   )*2.0
      V(2   )=V(4   )*2.0
      U(1   )=U(1   )+UDR
      U(2   )=U(2   )+UDI
C      WRITE(20,7777)'TERTIARY POTENTIALS'
C      WRITE(20,7777)'CENTRAL POTENTIAL'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)'SPIN ORBIT POTENTIAL'
C      WRITE(20,7778)(V(I),I=1,KT2)
C      WRITE(20,7777)'UD - DARWIN TERM'
C      WRITE(20,7778)(UD(I),I=1,KT2)
C 7777 FORMAT(A30)
C 7778 FORMAT(' ',1P10E12.4)
      ENDIF
C
      IF(IDIRAC.EQ.1) IDIRAC=0
      RETURN
C
 9000 FORMAT(10F8.4)
 9100 FORMAT(5E16.7)
 9508 FORMAT(18X,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4
     1          ,9H   2*S  =,F9.4,9H   FISW =,F9.4)
 9509 FORMAT(A18,9H   V RL =,F9.4,9H   R0RL =,F9.4,9H   A RL =,F9.4
     1          ,9H   R RL =,F9.4,9H   VSOR =,F9.4)
 9510 FORMAT(18X,9H   V IM =,F9.4,9H   R0IM =,F9.4,9H   A IM =,F9.4
     1          ,9H   R IM =,F9.4,9H   VSOI =,F9.4,9H   POWR =,F9.4)
 9512 FORMAT(18X,9H   BETA1=,F9.4,9H   LDFR1=,F9.4,9H   BETA2=,F9.4
     1          ,9H   LDFR2=,F9.4,9H   BETA3=,F9.4,9H   LDFR3=,F9.4)
 9515 FORMAT('0',20('*'),'  WARNING, Mixing of Dirac and non-Dirac'
     1,          ' potentials may be hazardous to your calculation'
     2,          20('*'))
      END