c$debug
c***********************************************************************
      SUBROUTINE ADWUCK4
c
c     Control program for first portion of DWUCK4 which reads in data
c     and calculates the form factors
c***********************************************************************
c
      parameter(ispc0 = 4010, ispc1 = 8000, ispc2 = 8000)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array0/space0(ispc0)
      Common/array1/space1(ispc1)
      Common/array2/space2(ispc2)
      DIMENSION U(800,2), V(800,2), UB(800), FF(10), VB(800)
      EQUIVALENCE (space0(   1), U ), (space0(1601), V)
     1           ,(space0(3201), UB), (space0(4001), FF) 
     2           ,(space1(   1), VB)
C
      IF(ICON(1).EQ.0) then
c
c       Use standard angle data
c
        II=(ANGLE(1)-ANGLE(2))/ANGLE(3)+1.0
        IF(II.LT.0) II=0
        ANGLE(1)=II
      else
C
C       READ CARD SET 2, Angle data
c
        READ (5,9001)ANGLE
      endif
      WRITE(6,9010)ANGLE
C
C
C     READ CARD SET 3, Angular momentum parameters
C
c     L      Maximum L for partial waves
c     Nltr   Number of angular momentum transfers
c     Ltrt   Orbital angular momentum transfer
c     Jtrt   Total   angular momentum transfer
c
      READ (5,9002)L,NLTR,(LTRT(I),I=1,NLTR),(JTRT(I),I=1,NLTR)
      WRITE(6,9008)L   ,(LTRT(I),I=1,NLTR)
      WRITE(6,9009)NLTR,(JTRT(I),I=1,NLTR)
      LC=0
      IF(L.LT.0) LC=IABS(L)
      L=IABS(L)
C
C
C     READ CARD SET 4
c
c     Drf    Radial increment
c     Rz     Lower radial cutoff
c     Rmaz   Upper Radial cutoff
c     Vce    Coulomb excitation scale factor for inelastic scattering reactions
c     Fnrng  Finite range parameter for stripping reactions
c     Amass  Mass for scaling radial increment, if zero, then scaling automatic
c
      READ (5,9001)DRF,RZ,RMAX,VCE,FNRNG,AMASS
      IF(ABS(RMAX).LT.ABS(RZ)) THEN
        temp=RZ
        RZ=RMAX
        RMAX=temp
      ENDIF
      WRITE(6,9011)DRF,RZ,RMAX,VCE,FNRNG
      IF(AMASS.NE.0.0) WRITE(6,9012)AMASS
C     default value for drf
      IF(DRF.EQ.0.0) DRF=0.1
      KZ=ABS(RZ  )/DRF+1.0E-08
      K =ABS(RMAX)/DRF+1.0E-08
      KC=2*((K+1)/2)
      LPLUS=L+1
      IND=0
      INCR=0
      IBF(5)=0
      IBF(6)=0
      WRITE(6,'(a)')' PARTICLE DATA '

C
C     READ IN DISTORTED WAVE INFORMATION
C
c     ibf(4) is flag to shut off spins if no spin orbit potential
c
      IBF(4)=0
c     Initial wave
      CALL FORMF(U(1,1),V(1,1),1,0,0,IVB)
      IF(IBF(4).NE.0) IBF(5)=1
      IBF(4)=0
c     Final wave
      CALL FORMF(U(1,2),V(1,2),2,0,0,IVB)
      IF(IBF(4).NE.0) IBF(6)=1
C
C     FORM FACTOR LOOP, up to 8 form factors may be used
C
      ILINE=0
      DO 35 II=1,NLTR
      IBF(4)=0
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      IF(ICON(3).NE.0.or.II.EQ.1) then
      DO 10 I=1,800
      UB(I)=0.0
   10 CONTINUE
      ivb=0
      do 11 I=1,10
      FF(I) = 0.0
   11 continue
c	IF
      IF(ICON(2).EQ.0.or.ICON(2).EQ.3) then
C
C     Here for collective or single particle form factors
C
      WRITE(6,'(a)')'0FORM FACTOR DATA '
      IF(RMAX.LT.0.0) K=KC
C
      CALL FORMF(UB,VB,3,LTRT(II),0,IVB)
      if(FK(3) .ne. 0.0) then
        FF(  7)=eta(3)*fk(3)
        FF(  8)=float(IVB)
        FF(  9)=float(kc)
      else
        FK(3) = sqrt(abs(FF(5)))
        IVB   = FF(8)
      endif
C
C     HERE FOR 2 NUCLEON FORM FACTORS
C
c	ELSEIF
      ELSEIF(ICON(2).EQ.2.OR.ICON(2).EQ.3) THEN
      CALL CATHEN(LTRT(II),JTRT(II),space2(   1),space2(1601),UB)
      ELSE
      write(6,'(a)')'0Warning, no form factor computed'
      endif
      endif
c
C     L.S FORM FACTOR KLUDGE
C     Move Spin-Orbit form factor to UB
      IF(IVB.EQ.1) THEN
        DO 14 M=1,800
        UB(M)=VB(M)
   14   CONTINUE
        FF(  8)=float(ivb)
      ENDIF
      ISTRT(II)=IS(3)
      IF(ICON(4).EQ.0.or.icon(4).eq.3) then
        IWRITE=2
        call ffprint(iwrite,k,alpha,ub(1),dr(3)
     1              ,ltrt(ii),istrt(ii),jtrt(ii),idat)
      endif
C
      IF(ICON(3).NE.0.or.II.EQ.1) then
C       Apply non-local correction to form factor
        CALL FNLOC(U,V,UB,VB,LTRT(II),Z,FNRNG,K)
      endif
      IF(VCE.NE.0.0) then
c     Store end point for Vcx (= Rmax*FF) for integration into complex R plane
      ivb = 4
        FF(  1)=3.0*VCE*CHSQ*Z(3)*ZA(3)
     1  *(RC(3)/(FLOAT(KC)*DR(3)))**LTRT(II)/FLOAT(2*LTRT(II)+1)
        FF(  2)=0.0
        FF(  3)=0.0
        FF(  4)=1.0
        FF(  5)=0.0
	FF(  6)=0.0
	FF(  7)=0.0
        FF(  8)=float(IVB)
        if(ltrt(ii).eq.0) ff(  1)=ff(  1)/3.0
      ENDIF
      IF(IVB.EQ.2.OR.IVB.EQ.4) THEN
C
C     SET ERROR TRAPS FOR UNBOUND STRIPPING AND COULOMB EXCITATION CASES
C
        ILINE=-1
C       RESTRICT MAX L TRANSFER IN UNBOUND STRIPPING - STORAGE LIMITED
        IF(LTRT(II).GT.7) THEN
          WRITE(6,9903)LTRT(II)
          IBF(3)=1
        ENDIF
        IF(FM(1).LT.FM(2)) THEN
C         SET ERROR FLAG
          IBF(3)=1
          WRITE(6,'(a)')'0ERROR * UNBOUND STATE FOR PICK UP NOT ALLOWED'
        ENDIF
C
        IF(IVB.EQ.2) THEN
          T1=(SQRT(FLOAT(LTRT(II)*(LTRT(II)+1))+ETA(3)**2)+ETA(3))/FK(3)
          IF(ABS(RMAX).LT.T1) WRITE(6,9904)T1,RMAX
        ENDIF
        L = FK(1)*FLOAT(KC)*DR(1) + 12.5
        IF(L.LT.LPLUS) WRITE(6,9905)
C
      ENDIF
C
C     STORE FF PARAMETERS FOR USE IN UNBND STRIP. AND Coulomb excitation
C
c     write UB and FF onto disk
      WRITE(4) UB,FF
      IS(3)=ISTRT(II)
      IF(IS(3).EQ.0) JTRT(II)=LTRT(II)+LTRT(II)
      IND=MAX0(LTRT(II)+1,IND)
      INCR=MAX0(INCR,JTRT(II))
C
      IF(ICON(4).EQ.3) then
        IWRITE=1
        call ffprint(iwrite,k,alpha,ub(1),dr(3)
     1              ,ltrt(ii),istrt(ii),jtrt(ii),idat)
      endif
C
C     CHECK ON SPINS AND STATISTICS
C
      LTR=LTRT(II)
      JTR=JTRT(II)
      IF(PHASEF(IS(1)+IS(2)+IS(3)).LT.0.0) GO TO 33
      IF(PHASEF(        JTR+IS(3)).LT.0.0) GO TO 33
      GO TO 35
   33 CONTINUE
C
C     SET ERROR FLAG
C
      IBF(3)=1
      WRITE(6,9906)IS,LTR,JTR
   35 CONTINUE
c
      IF(NLTR.GT.8) IBF(3)=1
C
C     CHECK ON COMPATIBILITY OF LMAX AND STORAGE
C
C			# Partial waves   # Radial matrix elements
      LPLUS=MIN0(LPLUS,400/(NS(1)+NS(2)),4000/(NS(1)*NS(2)*(IND+1)))
      L=LPLUS-1
      LPL2=LPLUS*2
      IF(RMAX.LT.0.0) K=KC
      IF(ILINE.NE.0) GO TO 40
      K=2*(K/2)
      L=FK(1)*float(K)*DR(1)+4.0
      L=MIN0(LPLUS-1,L)
      IF(LC.NE.0) L=MIN0(L,LC)
      LPLUS=L+1
      LPL2=LPLUS*2
   40 WRITE(6,9506)L,K,NLTR
      INC=(INCR+IS(1)+IS(2))/2+1
c                     # Plm's        # Fll's
      IBUF=MAX0(IBUF, INC*LPLUS+1, LPL2*IND+1)
      IBF(1)=IBUF
C
      IBF(7)=IS(1)
      IF(IBF(5).NE.0) GO TO 51
      IS(1)=0
      FS(1)=0.0
      NS(1)=1
   51 CONTINUE
      IBF(8)=IS(2)
      IF(IBF(6).NE.0) GO TO 52
      IS(2)=0
      FS(2)=0.0
      NS(2)=1
   52 CONTINUE
      IF(ICON(15).NE.0) THEN
C
C     PRINT OUT K(R)**2 FOR DISTORTED WAVES
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      DO 55 I=1,2
      DO 55 N=1,2
      IF(N.EQ.1) THEN
      WRITE(6,9057)I
      ELSE
      WRITE(6,9058)I
      ENDIF
      IND=0
      R=DR(I)
      DO 54 M=1,K,5
      MK=M+M-1
      MK4=MIN0(MK +9,K+K)
      IF(N.EQ.1) THEN
      WRITE(6,9052)R,(U(J,I),J=MK,MK4)
      ELSE
      WRITE(6,9052)R,(V(J,I),J=MK,MK4)
      ENDIF
      R=R+5.0*DR(I)
   54 CONTINUE
   55 CONTINUE
      ENDIF
   56 CONTINUE
      RETURN
C
 9001 FORMAT(10F8.4)
 9002 FORMAT(18I3)
 9008 FORMAT(18H0CARD SET 3 DATA  ,9H   LMAX =,I4,14H          LTR=,8I4)
 9009 FORMAT(18H                  ,9H   NLTR =,I4,14H        2*JTR=,8I4)
 9010 FORMAT(18H0ANGLE DATA       ,9H   THETN=,F9.4,9H   THET1=,F9.4
     1,           9H   DTHET=,F9.4,9H   A-ANG=,F9.4,9H   B-ANG=,F9.4)
 9011 FORMAT(18H0CARD SET 4 DATA  ,9H   DRF  =,F9.4,9H   RZ   =,F9.4
     1,           9H   RMAX =,F9.4,9H   VCE  =,F9.4,9H   FNRNG=,F9.4)
 9012 FORMAT(18X,9H   AMASS=,F9.4)
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9057 FORMAT(38H0Central    K(R)**2 FOR DISTORTED WAVE,I2)
 9058 FORMAT(38H0Spin orbit K(R)**2 FOR DISTORTED WAVE,I2)
 9500 FORMAT(1P5E16.7)
 9506 FORMAT(1H0,17X,9H   LMAX =,I4,14H        NSTEP=,I4
     1,                             14H        NLTR =,I4)
c 9803 FORMAT(1H ,16I2,4X,15A4,I2,2(1H/,I2.2),I4,2(1H.,I2.2))
 9903 FORMAT('0ERROR ** LTR=',I3,'  is too large' )
 9904 FORMAT('0*****  Warning  *****  CLASSICAL TURNING POINT FOR BOUND'
     1,'STATE =',F8.4,'  IS GREATER THAN RMAX =',F8.4)
 9905 FORMAT('0*****  Warning *****  LMAX too large for reliable '
     1,'calculation of unbound stripping or coulomb excitation cases.'
     2,'  Increase RMAX ')
 9906 FORMAT(1H0,28HSPIN STATISTICS NOT CORRECT ,7H 2*IS1=,I3,7H 2*IS2=
     1,I3,7H 2*IS3=,I3,7H   LTR=,I3,7H 2*JTR=,I3)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2),I4,2(1H.,I2))
      END

c***********************************************************************
      subroutine ffprint(ii,k,alpha,ub,drf,ltrt,istrt,jtrt,idat)
C
C     Prints form factor 
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ub(800),idat(6),alpha(15)
c
      if(ii.eq.1) then
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      endif
      R=DRF
      WRITE(6,9100)LTRT,ISTRT,JTRT
      WRITE(6,9101)
      DO 3006 M=1,K,5
      MK=M+M-1
      MK4=MIN0(MK+9,K+K)
      WRITE(6,9052)R ,(UB(N),N=MK,MK4)
      R=R+5.0*DRF
 3006 CONTINUE
      return
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9100 FORMAT(12H0FORM FACTOR,6X,9H    LTR =,I4,5X,9H   2*STR=,I4,5X
     1,                         9H   2*JTR=,I4)
 9101 FORMAT(55H0    R    RL,R        IM,R        RL,R+DR*1   IM,R+DR*1
     1,             48H   RL,R+DR*2   IM,R+DR*2   RL,R+DR*3   IM,R+DR*3
     2,             24H   RL,R+DR*4   IM,R+DR*4 )
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      end

c***********************************************************************
      SUBROUTINE FORMF(U,V,N,LAM,IK,IVB)
c
c     Processes kinematic input for waves and form factors and
c     calculates potentials or form factors
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
      COMMON/POTTER/DRX,AFACT(2),VFACT,SFACT,E,RM,G(4),ETAX,FKX,ETAKX
     1             ,RCX,HBARC2,ABETA(3),FLDF(3)
     2             ,NX,LAMX,KMXX,KX,IBX,LPLUSX,ICON4,NSPC,IDIRAC,ICHK
      DIMENSION U(810),V(800),X(4)
      EQUIVALENCE (X,FN)
C
C     READ IN POTENTIAL CARDS FOR CARD SETS 5,6,OR 7
C
      ETA6=60.
      KMAX=400
      ICON4=ICON(4)
      NSPC=N
      IVB=0
      KB=K
      READ (5,9000)E,FM(N),Z(N),FMA(N),ZA(N),RY,AC(N),PNLOC(N),FS(N),QCD
      E=E+QCD
      IF(N.EQ.2) THEN
        Q=E
        E=(ECM(1)+Q)*(FM(2)+FMA(2))/FMA(2)
      ENDIF
      IS(N)=FS(N)
      NS(N)=IS(N)+1
      IF(AMASS.EQ.0.0) AMASS=FMA(1)
      IF(IK.EQ.0) DR(N)=DRF*AMASS/FMA(N)
      KMXX=KMAX
      DRX=DR(N)
      AFACT(1)=FMA(N)**.333333333
      AFACT(2)=FM (N)**.333333333
      RC(N)=ABS(RY)*AFACT(1)
      IF(RY.LT.0.0) RC(N)=RC(N)+ABS(RY)*AFACT(2)
      DO 12 M=1,KMAX
      MK=M+M-1
      U(MK  )=0.0
      U(MK+1)=0.0
      V(MK  )=0.0
      V(MK+1)=0.0
   12 CONTINUE
      RM=0.0
      IF(E.EQ.0.0) GO TO 66
C
C        ICON(10).NE.0   GIVES RELATIVISTIC KINEMATICS
C
      IF(ICON(10).NE.0 .or. (fm(n).eq.0.0 .and. is(n).eq.2)) THEN
        IF(N.NE.2) GO TO 26
        IF(QCD.GT.0.0) GO TO 26
        E=E+(ECM(1)+Q)**2/(2.0*FMA(2)*AMU)
   26   CONTINUE
        FM1=FM(N)*AMU
        FM2=FMA(N)*AMU
        FMT=FM1+FM2
C     KLUDGE FAKE KE-LAB FOR BOUND STATES
        IF(N.GE.3) E=((E+FMT)**2-FMT**2)/(2.0*FM2)
        WLAB=E+FMT
        WCM =SQRT(2.0*E*FM2+FMT**2)
        GAMM=WLAB/WCM
        W1=GAMM*(FMT*FM1+FM2*E)/WLAB
        W2=GAMM* FM2
        ECM(N)=WCM-FMT
        IF(ICON(10).EQ.1) THEN
C       OLD KINEMATICS  SELECTED BY ICON(10) = 1
          FMU(N)=W1*W2/(W1+W2)
          VFACT=2.0*FMU(N)/HBARC**2
          SFACT=2.0*fm(n)*fma(n)/(fm(n)+fma(n))/hbarc**2
          EFACT=VFACT
          ELSE
C       NEW KINEMATICS  SELECTED BY ICON(10) > 1
          FMU(N)=W1
          VFACT= 2.0*W1 /HBARC**2
          SFACT= 2.0*FM1/HBARC**2
          EFACT=VFACT
        ENDIF
        FMU(N)=FMU(N)/AMU
        FK2(N)=(W1**2-FM1**2)/HBARC**2
      ELSE
C     NON RELATIVISTIC KINEMATICS
        FMU(N)=FM(N)*(FMA(N)/(FM(N)+FMA(N)))
        IF(N.LT.3) ECM(N)=E*(FMA(N)/(FM(N)+FMA(N)))
        IF(N.GE.3) ECM(N)=E
        VFACT=2.0*FMU(N)*AMU/HBARC**2
        SFACT=VFACT
        EFACT=VFACT
        FK2(N)=SFACT*ECM(N)
      ENDIF
      FK(N)=SQRT(ABS(FK2(N)))
      ETAK=CHSQ*Z(N)*ZA(N)*EFACT
      ETA(N)=ETAK*0.5/FK(N)
      HBARC2=HBARC**2
C
C     ADD COULOMB AND KINETIC ENERGIES TO U
C
      RCX=RC(N)
      IF(RCX.EQ.0.0) RCX=DR(N)
      R=0.0
      FCOU=0.5*ETAK/RCX
      DO 42 M=1,kmax
      MK=M+M-1
      R=R+DR(N)
      IF(R.GT.RCX) GO TO 40
      FC=FCOU*(3.0-(R/RCX)**2)
      GO TO 41
   40 CONTINUE
      FC=ETAK/R
   41 CONTINUE
      IF(N.NE.3) U(MK  )=U(MK  )+FK2(N)-FC
      IF(N.EQ.3) U(MK+1)=U(MK+1)+FK2(N)-FC
   42 CONTINUE
      GO TO 67
   66 CONTINUE
      FK(N)=0.0
      ETA(N)=0.0
      FK2(N)=0.0
      ECM(N)=0.0
      FMU(N)=FM(N)*FMA(N)/(FM(N)+FMA(N))
      VFACT=2.0*FMU(N)*AMU/HBARC**2
      SFACT=VFACT
      ETAK=CHSQ*Z(N)*ZA(N)*VFACT
   67 CONTINUE
      IF(N.LT.3) GO TO 68
      IF(ICON(4).EQ.2) GO TO 69
   68 CONTINUE
      IF(N.GT.2) THEN
        Q=ECM(N)
      ELSE
        Q=ECM(N)-ECM(1)
      ENDIF
      WRITE(6,9010)N
      WRITE(6,9503)E,RY,AC(N),FS(N)
      WRITE(6,9504)FM(N),FMA(N),Q
      WRITE(6,9505)Z(N),ZA(N),PNLOC(N)
      WRITE(6,9500)
      RHO=FK(N)*RC(N)
      WRITE(6,9506)ECM(N),RC(N),RHO
      WRITE(6,9507)FK(N),ETA(N),DR(N)
      WRITE(6,9008)
   69 CONTINUE
      FS(N)=FS(N)/2.0
C     Set up variables for potential routine
      IBX=0
      ETAX=ETA(N)
      FKX=FK(N)
      ETAKX=ETAK
      RCX=RC(N)
      LAMX=LAM
      NX=N
      LPLUSX=LPLUS
      KX=K
      ICHK=0
      IDIRAC=0
C
      CALL POTS(U,V)
C
      DR(N)=DRX
      LPLUS=LPLUSX
      IBF(4)=IBX
C     Set nonlocality for Dirac-Darwin term
      IF(N.LE.2.AND.IDIRAC.NE.0.AND.PNLOC(N).EQ.0.0) PNLOC(N)=-1.0
      K=MIN0(MAX0(K,KX),KMAX)
      IF(N.LE.2) GO TO 3000
      IF(E.EQ.0.0) THEN
      IF(IBF(4).NE.0) u(808) = 1
      IF(IVB.EQ.4) THEN
        U(801)=U(2*KC-1)*FLOAT(KC)*DR(3)
        U(802)=0.0
        U(803)=0.0
        U(804)=1.0
      ENDIF
      GO TO 3000
      ENDIF
C
C     SINGLE PARTICLE ORBITAL
C
C
C     READ IN QUANTUM NUMBERS FOR SINGLE PARTICLE ORBITAL
C
 2000 CONTINUE
C     Set flags for unbound stripping case
      IF(E.GT.0.0) THEN
        IVB=2
        K=KC
      ENDIF
      READ(5,9000)G,VTRIAL,FISW
      FN =G(1)
      FL =G(2)
      FJ2=G(3)
      FSS=G(4)
      ISW=FISW
      WRITE(6,9500)
      FJ0=FJ2/2.0
      FS0=FSS/2.0
      IF(VTRIAL.EQ.0.0) VTRIAL=ETA6
      WRITE(6,9508)G,FISW
      FACT=(FJ0**2+FJ0-FL**2-FL-FS0**2-FS0)*0.5
c
      DO 2028 M=1,kmax
      MK=M+M-1
      V(MK  )=U(MK  )+V(MK  )*FACT
      V(MK+1)=U(MK+1)+V(MK+1)*FACT
 2028 CONTINUE
      WRITE(6,9500)
c
      CALL BIND(V,  U     ,DR(3),RM,FN,FL,K,FK(3),ETA(3),VTRIAL,ECM(3)
     1,FK2(3),ISW,IBF(3),U(801))
C
      IBF(2)=RM/DR(3)
      DO 2050 M=1,K
      MK=M+M-1
      V(M)=VTRIAL*V(MK  )+V(MK+1)
 2050 CONTINUE
      Anorm=1.0
      FACT=PNLOC(3)**2/8.0
      IF(FACT.ne.0.0) then
C
C     NON-LOCAL CORRECTION FOR SINGLE PARTICLE FUNCTION
C
      SUM=0.0
      R=0.0
      DO 2075 M=1,K
      MK=M+M-1
      R=R+DR(3)
      U(M)=U(M)*EXP(FACT*(FK2(3)-V(M)))
      SUM=SUM+(U(M)*R)**2
 2075 CONTINUE
C
C     DO NOT RENORMALIZE FOR POSITIVE ENERGY
C
      IF(FK2(3).LT.0.0) then
        Anorm=1.0/SQRT(SUM*DR(3))
      ELSE
        Anorm=1.0
      ENDIF
      ENDIF
      DO 2100 M=kmax,1,-1
      MK=M+M-1
      IF(M.GT.K) then
        V(M   )=0.0
        U(MK  )=0.0
      ELSE
        U(MK  )=U(M)*Anorm
      ENDIF
      U(MK+1)=0.0
 2100 CONTINUE
c
 3000 CONTINUE
c     Store quantum numbers for transfer back
      DO 3020 M=1,4
      X(M)=G(M)
 3020 CONTINUE
      RETURN
 9000 FORMAT(10F8.4)
 9008 FORMAT(21H0POTENTIAL PARAMETERS )
 9010 FORMAT( 9H0PARTICLE,I2,115(1H*))
 9500 FORMAT(1H ,3A6,5(3X,A6,F9.4))
 9503 FORMAT(18H INPUT DATA       ,9H   ELAB =,F9.4,9H   RC0  =,F9.4
     1           ,9H   AC   =,F9.4,9H   2*STR=,F9.4)
 9504 FORMAT(18X ,9H   MASSP=,F9.4,9H   MASST=,F9.4,9H   Q    =,F9.4)
 9505 FORMAT(18X ,9H   ZP   =,F9.4,9H   ZT   =,F9.4,9H   PNLOC=,F9.4)
 9506 FORMAT(18H DERIVED DATA     ,9H   ECM  =,F9.4,9H   RC   =,F9.4
     1           ,9H   RHO  =,F9.4)
 9507 FORMAT(18X ,9H   K    =,F9.4,9H   ETA  =,F9.4,9H   DR   =,F9.4)
 9508 FORMAT(18X ,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4
     1           ,9H   2*S  =,F9.4,9H   FISW =,F9.4)
 9511 FORMAT(18X ,9H   DAMP =,F9.4)
      END

c***********************************************************************
      SUBROUTINE FNLOC(U,V,W,VB,LTR,BKIN,FNRNG,KT)
c
c     Calculates non-local and finite range correction for distorted waves
c
c     U     = Optical     potentials for distorted waves
c     V     = Spin orbit  potentials for distorted waves
c     W     = Form factor
c     VB    = Bound State potential
c     LTR   = Angular momentum transfer ( of bound state)
c     BKIN  = Kinematic quatities stored in blank common
c     FNRNG = Finite range parameter
c     KT    = Number of points
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(800,2),V(800,2),W(800),VB(800),BKIN(3,14)
     1,RT(8,7),FM(3),FMA(3),RC(3),PNLOC(3),FMU(3),FK(3),FK2(3)
     1  ,ETA(3),DR(3)
      DATA  RT/  0.0   ,  4.4934,  7.7253, 10.9041, 14.0662, 17.2208,
     1 20.3713, 23.5195,  2.0816,  5.9404,  9.2058, 12.4044, 15.5792,
     2 18.7426, 21.8997, 25.0528,  3.3421,  7.2899, 10.6139, 13.8461,
     3 17.0429, 20.2219, 23.3905, 26.5526,  4.5141,  8.5838, 11.9727,
     4 15.2445, 18.4681, 21.6666, 24.8501, 28.0239,  5.6467,  9.8404,
     5 13.2956, 16.6093, 19.8624, 23.0828, 26.2833, 29.4706,  6.7565,
     6 11.0702, 14.5906, 17.9472, 21.2311, 24.4748, 27.6937, 30.8960,
     7  7.8511, 12.2793, 15.8632, 19.2627, 22.5781, 25.8461, 29.0843,
     8 32.3025/
C
C     FINITE RANGE CORRECTION
C
C
C     POSITIVE FR PARAMETER CALCULATES A HULTHEN TYPE FR CORRECTION
C     NEGATIVE FR PARAMETER CALCULATES A EXPONENTIAL TYPE FR CORRECTION
C
      K=KT
      DO 20 I=1,3
      FM(I)   =BKIN(I,3)
      FMA(I)  =BKIN(I,4)
      RC(I)   =BKIN(I,5)
      PNLOC(I)=BKIN(I,7)
      FK(I)   =BKIN(I,10)
      FK2(I)  =BKIN(I,11)
      ETA(I)  =BKIN(I,12)
      DR(I)   =BKIN(I,13)
      FMU(I)  =BKIN(I,14)
   20 CONTINUE
      IF(FNRNG.EQ.0.0) GO TO 398
      IF(FM(1)-FM(2))143,398,142
  142 I1=1
      I2=2
      GO TO 144
  143 I1=2
      I2=1
C
C     I1 IS FOR LARGER MASS , I2 IS FOR SMALLER MASS
C
  144 CONTINUE
      FM1=FMU(I1)
      FM2=FMU(I2)
      TEMP=FM(I2)*FNRNG**2/FM(I1)
      FACT=TEMP*FMA(I1)/FMA(I2)
      FMX=FM(I1)-FM(I2)
      TEMP=TEMP*FMX
      IF(FMX.LT.1.9) GO TO 350
C
C     FORM FACTOR PART BY ROST-KUNZ METHOD
C     THIS OPTION CALCULATES AN EXPONENTIAL TYPE FR CORRECTION
C     FOR MULTI -NUCLEON TRANSFER REACTIONS
C
      DR2=DR(3)**2
      FLF=LTR*(LTR+1)
      R3=DR(3)*FLOAT(K+1)
      D3=0.0
      A3=0.0
      KM=K+K+1
      R2=R3-DR(3)
      D2=-FLF/R3
      A2=W (KM-2)*R2
      DO 220 M=1,K
      KM=KM-2
      R1=R2-DR(3)
      A1=W (KM-2)*R1
      D1=-FLF/R2**2
      IF(A3.EQ.0.0) GO TO 215
      D1=D1+(A3-2.0*A2+A1)/(A2*DR2)
      W (KM+2)=W (KM+2)*EXP(D2*FACT)
      DX2=(D3*A3-2.0*D2*A2+D1*A1)/(R2*DR2)-A2*D2**2/R2
      W (KM+2)=W (KM+2)+0.5*DX2*FACT**2
  215 CONTINUE
      IF(ABS(W (KM-2)).LT.ABS(W (KM  ))) GO TO 225
      R3=R2
      R2=R1
      D3=D2
      D2=D1
      A3=A2
      A2=A1
  220 CONTINUE
  225 CONTINUE
      KM=KM/2+1
      DO 227 M=1,KM
      MK=M+M-1
      VB(M)=W(MK)
      W(MK)=0.0
  227 CONTINUE
      RKM=FLOAT(KM)*DR(3)
      LL=MIN0(LTR,6)
      FLF=LL*(LL+1)
      DR12=DR(3)**2/12.0
      DO 250 N=1,8
      FKX=(RT(N,LL+1)/RKM)**2
      A1=0.0
      D1=0.0
      R=DR(3)
      A2=DR(3)
      D2=(1.0+DR12*(FKX-FLF/R**2))*A2
      FNORM=0.0
      AN=0.0
      DO 230 M=1,KM
      AN=AN+VB(M)*A2*R
      FNORM=FNORM+A2**2
      VB(M+400)=A2/R
      R=R+DR(3)
      D3=1.0+DR12*(FKX-FLF/R**2)
      A3=(12.0*A2-10.0*D2-D1)/D3
      D1=D2
      D2=D3*A3
      A1=A2
      A2=A3
  230 CONTINUE
      AN=AN-0.5*VB(KM)*A1*(R-DR(3))
      FNORM=FNORM-0.5*A1**2
      FNORM=AN*EXP(-FACT*FKX)/FNORM
      DO 240 M=1,KM
      MK=M+M-1
      W (MK)=W (MK)+VB(M+400)*FNORM
  240 CONTINUE
  250 CONTINUE
  340 CONTINUE
C
C     DISTORTED WAVE PART, GAUSSIAN FACTOR
C
      DO 345 M=1,K
      MK=M+M-1
      CTEMP1=TEMP*(U(MK  ,I2)/FM2-U(MK  ,I1)/FM1)
      CTEMP2=TEMP*(U(MK+1,I2)/FM2-U(MK+1,I1)/FM1)
      FACT=EXP(-CTEMP1)
      UT1= FACT*COS(CTEMP2)
      UT2=-FACT*SIN(CTEMP2)
      UF1   =UT1*W(MK  )-UT2*W(MK+1)
      UF2   =UT1*W(MK+1)+UT2*W(MK  )
      W(MK  )=UF1
      W(MK+1)=UF2
  345 CONTINUE
  348 CONTINUE
      GO TO 398
  350 CONTINUE
      IF(FNRNG.GT.0.0) GO TO 360
C
C     SINGLE PARTICLE PART, GAUSSIAN FACTOR
C
      DO 355 M=1,K
      MK=M+M-1
      FACT=EXP(-TEMP*VB(M)/FMU(3))
      W(MK  )=W(MK  )*FACT
      W(MK+1)=W(MK+1)*FACT
  355 CONTINUE
      GO TO 340
  360 CONTINUE
C
C     HULTHEN TYPE
C
      DO 365 M=1,K
      MK=M+M-1
      CTEMP1=VB(M)/FMU(3)+U(MK  ,I2)/FM2-U(MK  ,I1)/FM1
      CTEMP2=             U(MK+1,I2)/FM2-U(MK+1,I1)/FM1
      CTEMP1=1.0+TEMP*CTEMP1
      CTEMP2=    TEMP*CTEMP2
      DET=CTEMP1**2+CTEMP2**2
      UT1=(CTEMP1*W(MK  )+CTEMP2*W(MK+1))/DET
      UT2=(CTEMP1*W(MK+1)-CTEMP2*W(MK  ))/DET
      W(MK  )=UT1
      W(MK+1)=UT2
  365 CONTINUE
  398 CONTINUE
C
C     NON LOCALITY CORRECTION FACTOR FOR DISTORTED WAVES
C
      DO 410 I=1,2
      CALL FNLOC5(U(1,I),V(1,I),W(1),PNLOC(I),FK2(I),ETA(I)
     1,           FK(I),RC(I),DR(I),K)
  410 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE CATHEN(LOPT,JOPT,UB,VB,SI)
c
c     Two particle form factor subroutine, calculates microscopic
c     inelastic or two-nucleon transfer form factors
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
      DIMENSION UB(800,2),VB(801),SI(800),C(10),QNUM(4,2),G(4),IQN(3)
     1         ,JI(2),JJ(2),FJT(2),CON(8)
      EQUIVALENCE (G,FN),(C(1),CNTROL),(C(2),QCODE),(C(3),FMUV)
     1,(C(4),OPT)
      DATA SQR4PI,SQRTEN,PI/3.54490780,3.162277660,3.141592/
C
      IQFLG=0
      KMFLG=0
      IQN(1)=LOPT
      IQN(3)=JOPT
      KT=399
      IBF(2)=1
      FK2(3)=0.0
  100 CONTINUE
      READ(5,9001)C
      CNTRL=ABS(CNTROL)
      IF(CNTRL.EQ.0.0) GO TO 3010
      IF(QCODE)2020,3000,1100
 1100 CONTINUE
      icode=abs(qcode)
C
C     LOOP FOR TWO ORBITALS
C
      ML=1
      DO 2015 I=1,2
      FJT(1)=C(5)
      FJT(2)=C(I+5)
      IF(I.EQ.2.AND.CNTRL.EQ.1.0) GO TO 2005
      IK=0
      IF(icode.EQ.7.OR.icode.EQ.8) ik=i-1
      G(1)=0.0
      G(2)=0.0
      G(3)=0.0
      G(4)=0.0
C
      CALL FORMF(UB(1,2),VB,3,LOPT,IK,IVB)
      IVB=0
C
      DRX=DR(3)
      IF(KMFLG.EQ.0) KM=IBF(2)
      KMFLG=1
      KT=MIN0(KT,K)
 2005 CONTINUE
      IF(FJT(1).EQ.0.0) FJT(2)=FJ2
C
C     PRINT OUT SINGLE PARTICLE INFORMATION
C
      if(icode .ne. 5) then
        WRITE(6,9501)
        WRITE(6,9100)FJ2,FJT
      endif
      IF(ICON(4).ne.2 .and. (CNTRL.EQ.2.0 .or. I.EQ.1)) then
        R=DRX
        K2=K+K
        DO 2006 M=1,K2,20
        MK4=MIN0(M+19,K2)
        WRITE(6,9052)R,(UB(N,2),N=M,MK4,2)
        R=R+10.0*DRX
 2006   CONTINUE
      endif
      DO 2008 M=1,400
      MK=M+M-1
      UB(ML,1)=UB(MK,2)
      ML=ML+1
 2008 CONTINUE
      DO 2009 M=1,4
      QNUM(M,I)=G(M)
 2009 CONTINUE
      JI(I)=FJT(1)
      JJ(I)=FJT(2)
 2015 CONTINUE
C
 2020 ICODE=ABS(QCODE)
      KT=MIN0(KT,K)
      LTR1=IQN(1)+1
      IQN(2)=IS(3)
      IF(IS(3).EQ.0) IQN(3)=IQN(1)+IQN(1)
      WRITE(6,9051)QCODE,FMUV,OPT
C     ENTER TIME REVERSAL PHASE
      LVR=QNUM(2,1)+QNUM(2,2)
      OPT=OPT*PHASEF((LVR+IQN(1))/2)
      IFLAG=1
      IF(ICODE.LE.10) GO TO 2025
      IFLAG=2
      ICODE=ICODE-10
 2025 CONTINUE
C
      GO TO (2100,2100,2300,2990,2500,2100,2700,2700),ICODE
C
C     ICODE=1  YUKAWA POTENTIAL
C     ICODE=2  COULOMB POTENTIAL
C     ICODE=3  OPEP TENSOR POTENTIAL
C     ICODE=5  TWO PARTICLE TRANSFER
C     ICODE=6  ZERO RANGE KNOCK OUT
C     ICODE=7  NO RECOIL FORM FACTOR
C
 2100 LVR=0
      TOPT=1.0
      GO TO 2350
 2300 ICODE=1
      LVR=2
      TOPT=-SQRTEN
 2350 MNNL=IABS(LTR1-LVR-1)+1
      MXXL=LTR1+LVR
      LTR2=LTR1+LTR1-2
      LVR2=LVR+LVR
      LLX=QNUM(2,1)+QNUM(2,1)
      LPX=QNUM(2,2)+QNUM(2,2)
      JLX=QNUM(3,1)
      JPX=QNUM(3,2)
      IS1=QNUM(4,1)
      IS2=QNUM(4,2)
      OPT=OPT*PHASEF(LPX/2)
      DRX=DR(1)
      DO 2390 LAM =MNNL,MXXL,2
      LAM2 =LAM +LAM -2
      VOPT=TOPT*OPT*SQR4PI*SQRT(FLOAT(LAM2 +1)*FLOAT(IQN(2)+1))
     1 *VCC(LAM2 ,LVR2,LTR2 ,0,0)
     2 *RACAH(IQN(2),IQN(2),LAM2 ,LTR2 ,LVR2,IQN(3))
      RME=0.0
      IF(IQN(2).EQ.0) RME=1.0
      IF(IQN(2).EQ.2) RME=SQRT(FLOAT(IS(1)*(IS(1)+2))
     1*FLOAT(IS1*(IS1+2)))
      RME=RME*SQRT(FLOAT(IS1+1))*PHASEF((IS(1)-IS(2))/2)
     1 *SQRT(FLOAT((LLX+1)*(LAM2+1)))*VCC(LLX,LAM2,LPX,0,0)
     2 *SQRT(FLOAT(JLX+1)*FLOAT(JPX+1)*FLOAT(IQN(3)+1))
       rme =rme*WINEJ(LPX,LLX,LAM2 ,IS2,IS1,IQN(2),JPX,JLX,IQN(3))
     4 *PHASEF((JI(1)+IQN(3)-JLX-JJ(2))/2)
     5 *RACAH(JPX,JJ(2),JLX,JJ(1),JI(1),IQN(3))
     6 *SQRT(FLOAT(JJ(1)+1)*FLOAT(JJ(2)+1))
      VOPT=VOPT*RME*SQRT(FLOAT(LTR2+1)/FLOAT((JJ(1)+1)*(IQN(3)+1)))
     1 *(DR(3)/DR(1))**3
      SL=0.0
      IF(VOPT.EQ.0.0) GO TO 2390
      IF(ICODE.EQ.6 ) GO TO 2365
C
      KMAX=400
      CALL SLATR (KT,KMAX,DRX,VB(  1),LTR1,FMUV,ICODE)
      CALL SLATR (KT,KMAX,DRX,VB(801),LAM ,FMUV,ICODE)
C
C
      CALL RADIN
     1(KT,KMAX,DRX,FMUV,VB,UB(1,1),UB(401,1),SL,VOPT,SI(IFLAG),KM,SK)
C
      GO TO 2380
 2365 CONTINUE
      SCALE=VOPT/SQR4PI**2
      R=0.0
      MK=IFLAG
      DO 2370 M=1,KT
      R=R+DRX
      TEMP=UB(M    ,1)*UB(M+400,1)*SCALE
      SI(MK  )=SI(MK  )+TEMP
      SL      =SL      +TEMP*R**2
      IF(M.EQ.KM) SK=TEMP
      MK=MK+2
 2370 CONTINUE
      SL=SL*DRX
 2380 CONTINUE
 2385 I=LAM -1
      WRITE(6,9102)SL,IQN(1),I,RME
      WRITE(6,9002)KM,SK
 2390 CONTINUE
      GO TO 2990
C
C     HERE FOR TWO NUCLEON TRANSFER
C
 2500 CONTINUE
c
c     c1 = R1 scale
c     c2 = R2 scale
c     c3 = r1 scale
c     c4 = r2 scale
c     c5 = r0, integration scale length for relative coordinate
c     c6 = Pauli flag
c     c7 = order  of gaussian integration
c     c8 = number of integration points
c
      IF(FMUV.EQ.0.0) FMUV=1.7
      CON(1)= 1.0
      CON(2)= 1.0
      CON(3)= 0.5
      CON(4)=-0.5
      CON(5)= 2.0*fmuv
      con(6)= 0.0
      con(7)= 0.5
      con(8)= 0.0
      OPT=OPT*16./PI
      CALL DSTRIP
     1(IQN,DRX,KT,UB(1,1),UB(401,1),SI(IFLAG),QNUM,OPT,KM,SL,CON)
      WRITE(6,9002)KM,SL
      IQFLG=1
      GO TO 2990
 2700 CONTINUE
      T1= MIN (FM(1),FM(2))
      T2= MAX (FM(1),FM(2))
      CON(1)= 1.0
      CON(2)= 0.0
      CON(3)= T1/T2
      CON(4)= 1.0
      CON(5)= 1.0
      CON(7)= 0.0
      CON(8)= 0.0
      IF(ICODE.EQ.8) GO TO 2706
      T3=HBARC**2/(2.0*AMU*FMU(3))
      R=0.0
      DO 2705 M=1,KT
      R=R+DRX
      UB(M+400,1)=UB(M+400,1)*T3*  (FK2(3)-VB(M))
 2705 CONTINUE
 2706 CONTINUE
      M=(KT+KM)/2+400
      T1= LOG(UB(M-1,1)*FLOAT(M-401)/(UB(M  ,1)*FLOAT(M-400)))
      FMUV = DRX/T1
      CON(5)=fmuv
      DO 2710 M=1,KT
      UB(M+400,1)=UB(M+400,1)*EXP(T1*FLOAT(M))*DRX*FLOAT(M)/100.
 2710 CONTINUE
      OPT=OPT*FMUV*FMUV*SQR4PI/2.0
      CALL DSTRIP
     1(IQN,DRX,KT,UB(1,1),UB(401,1),SI(IFLAG),QNUM,OPT,KM,SL,CON)
      WRITE(6,9002)KM,SL
      IQFLG=1
      GO TO 2990
 2990 CONTINUE
 3000 CONTINUE
      IF(CNTROL.GT.0.0) GO TO 100
      IF(IQFLG .EQ.1) KT=FLOAT(KT)-1.50/DRX
      K=KT
 3010 RETURN
 9001 FORMAT(10F8.4)
 9002 FORMAT(15H FORM FACTOR,M=,I3,1PE18.6)
 9051 FORMAT(18H0PARAMETERS       ,9H   QCODE=,F9.4,9H   RANGE=,F9.4
     1 ,9H   VZERO=,F9.4)
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9100 FORMAT(27H0SINGLE PARTICLE FUNCTIONS ,8HCOUPLING ,6H 2*J1=,F3.0
     1 ,9H    2*J2=,F3.0,9H    2*JI=,F3.0)
 9102 FORMAT(19H VOLUME INTEGRAL = ,F10.4,8H    LTR=,I3,8H    LAM=,I3
     1,8H    RME=  ,F8.4)
 9501 FORMAT(1H ,18HFORM FACTOR DATA  )
      END
