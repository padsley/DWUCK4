
c$debug
c***********************************************************************
      SUBROUTINE RADINT(F,FLL,UB,FF,LTR)
c
c     Subroutine for computing radial integrals
c
c     f     Radial wave function storage
c     fll   Radial integral storage
c     ltr   Angular momentum transfer
c***********************************************************************
c
c     Ymax   maximum X = ikr, off real axis
c     Nx     number of points off real axis, must be a multiple of 5
C
      parameter(ispc0 = 2005, ispc1 = 4000, ispc2 = 4000)
      IMPLICIT REAL*8(A-H,O-Z)
      double complex space0, space1, space2
      parameter (kmax = 400, lmax = 400, nx=200, Ymax=10.0)
      double complex FLL(*), F(*), H(lmax), C(lmax), D(lmax)
     1  , cfk1, cfk2, cfk3, FP(2*nx+1), FA(2*nx+1)
     1  , UB, CTEMP, VTEMP, FFTEMP, buffr
     2  ,  T1, CN, CP, CD, SUM1, SUM3, Smat, HN
c
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
c
      Common/array0/space0(ispc0)
      Common/array1/space1(ispc1)
      Common/array2/space2(ispc2)
c
      DIMENSION UB(kmax),ff(10),buffr(1600)
     1,         DRX(2), LPTBL(lmax/2), S(lmax), WT(5)
      EQUIVALENCE (space1(   1), buffr), (space0(   1), D)
     2           ,(space0( 401), H), (space0(1201), C)
C
      IF(RC(3).EQ.0.0) RC(3)=DR(3)
c     Coulomb excitation factor
      VFACT=3.0*VCE*CHSQ*Z(3)*ZA(3)*RC(3)**LTR/FLOAT(LTR+LTR+1)
      if(ltr.eq.0) vfact=vfact/3.0
      DRX(1)=2.0*DR(3)/3.0
      DRX(2)=2.0*DRX(1)
C     JACOBIAN FACTOR
      ANORM= MAX (DR(1),DR(2))**2/(DR(3)* MIN (DR(1),DR(2)))
      JR = NS(1)
      JS = NS(2)
      NP=LPLUS*(JR+JS)
      INCR=LPLUS*(LTR+1)
c
C     Unbound state and coulomb excitation end point data
C     ff( 1) = D  REAL  (exp(2i*delta) - 1)/(2i)
C     ff( 2) = D  IMAG  (exp(2i*delta) - 1)/(2i)
C     ff( 3) = G  irregular solution at R = Rmax
C     ff( 4) = F  regular   solution at R = Rmax
C     ff( 5) = K**2 real
C     ff( 6) = K**2 imag
C     ff( 7) = Eta * K
C     ff( 8) = IVB switch
C		= 2, unbound state routine
C		= 4, coulomb excitation
C     ff( 9) = Number of last point of form factor
C     ff(10) = ISW parameter in BIND routine
C
      IVB   =ff( 8)
      isw   =ff(10)
c
      if(ivb.eq.2) then
c       Unbound state
        FFTEMP = UB(KC)
      elseif(ivb.eq.4) then
c       Coulomb excitation
        fftemp = cmplx(ff(1), ff(2))/(dr(3)*float(k))
      endif
c     Spin orbit factor
      FSFACT=FS(1)*(FS(1)+1.)+FS(2)*(FS(2)+1.)
c     clear radial integral storage
      INC=INCR*NS(1)*NS(2)
      DO 3 II=1,INC
      FLL(II)=0.0
    3 CONTINUE
c
      R3 = 0.0
      MX = 1
      DO 100 M=1,K
      MX = 3-MX
      R3 = R3+DRF
C
C     READ IN DISTORTED WAVE RADIAL FUNCTIONS
C
      READ (4)(F(J),J=1,NP)
      VTEMP = UB(M)*DRX(MX)
      IF(VFACT.NE.0.0) then
C
C     COULOMB EXCITATION ADDITION TO FORM FACTOR
C
      IF(R3 .GT .RC(3)) then 
      VTEMP = VTEMP + DRX(MX)*VFACT/R3**(LTR+1)
      endif
      endif
C
C**** READ IN NORMALIZATION CONSTANTS at last point
C
      IF(M .EQ. K) then
        READ (4)(C(I),I=1,NP)
c       adjust last point in integration rule
        VTEMP = VTEMP*0.5
      endif
c
      INC=0
      IY=0
c
c     Loop over initial state spin states
c
      FS1=-FS(1)
      DO 96 I= 1,JR
      IZ=NS(1)*LPLUS
c
c     Loop over final   state spin states
c
      FS2=-FS(2)
      DO 95 J=1,JS
c
c     Final  state L values
c
      DO 90 LL=1,LPLUS
      LK=LL
      IZ1=LL+IZ+(J-1)*LPLUS
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      CTEMP = VTEMP*F(IZ1)
c
c     Initial state L values
c
      DO 80 LP=LP1,LP2,2
      IY1=LP+(I-1)*LPLUS
      IND=LK+INC
      IF(M .LT. KZ) GO TO 79
      FLL(IND)=FLL(IND) + CTEMP*F(IY1)
      IF(M.EQ.K) THEN
C
C     NORMALIZE RADIAL INTEGRALS
C
      FLL(IND) = C(IZ1)*FLL(IND)*C(IY1)*ANORM
C     L*S FORM FACTOR,  use symmetrized form
      IF(IVB.EQ.1) THEN
        FJ1=FLOAT(LL)+FS1
        FJ2=FLOAT(LP)+FS2
        TEMP1 = (FJ1*(FJ1-1.)+FJ2*(FJ2-1.)-FLOAT(LL*(LL-1))
     1        -FLOAT(LP*(LP-1))-FSFACT)/4.
        FLL(IND)=FLL(IND)*TEMP1
      ENDIF
      ENDIF
   79 continue
      LK=LK+LPLUS
   80 CONTINUE
   90 CONTINUE
      INC=INC+INCR
      FS2=FS2+1.0
   95 CONTINUE
      IY=IY+LPLUS
      FS1=FS1+1.0
   96 CONTINUE
  100 CONTINUE
c     special flag to skip this section
      IF(ICON(20).NE.0) GO TO 1100
      IF(IVB.EQ.2 .OR .IVB.EQ.4) THEN
C
C     UNBOUND STATE SECTION
C
      cfk1 =      fk2(1)
      cfk2 =      fk2(2)
      cfk3 =cmplx(ff (5),ff(6))
      CN   =cmplx(FF(1), FF(2))
      HN   =cmplx(FF(3), FF(4))
      etak  =ff(7)
      FK(3) = sqrt(abs(ff (5)))
c
C     read in scatt. amp. and asymptotic fcts - h(ODD) = g, h(EVEN) = f
c     and coulomb phases
c
      READ (4)(D(I),I=1,NP),(H(I),I=1,LPL2),(S(I),I=1,LPL2)
      R1=DR(1)*FLOAT(KC)
      R2=DR(2)*FLOAT(KC)
      R3=DR(3)*FLOAT(KC)
      DX=Ymax/FLOAT(NX)
      EZ=EXP(-DX)
      DX1=DX/(FK(1)-FK(2)*FMA(1)/FMA(2)-FK(3))
      IF(DX1.LE.0.0) THEN
        WRITE(6,'(a)')'0KINEMATICS INCONSISTENT - EXIT ON ERROR'
        STOP 'Kinematics error in Routine UNBIND'
      ENDIF
C
      NI=DX1/DR(1)/2.0+0.5
      NI=MAX0(NI,1)
      DX1=DX1/FLOAT(NI)
      DX2=DX1*FMA(1)/FMA(2)
C     Calculate form factor extension into complex (Rmax,y) plane
      IFLG=1
C     -----------------------------------------------------
      CALL UNBIND(LTR,NX,NI,FA,CN,HN,etak,cfk3,DX1,R1,IFLG)
C     -----------------------------------------------------
C     ADJUST FA FOR FINITE RANGE, NON LOCALITY, AND PHASE
C     FFTEMP IS FORM FACTOR AT RMAX
C     FA(R+i0) IS UNMODIFIED FORM FACTOR FROM BIND AT RMAX
C     T1 = FFTEMP/FA(R+i0)
      T1 = FFTEMP/FA(nx+1)
      DX=DX1*FLOAT(NI)
      NX1=NX+1
C     Multiply FA BY R/(R+iy )
      DO 120 M=1,NX1 
      M1=M
      M2=2*NX+2 - M
      X  = dx*float(NX1-M) 
      FA(M1) = FA(M1)*(R1*T1)/cmplx(R1, X)
      FA(M2) = conjg(FA(M1))
  120 CONTINUE
C
      WT(1)=190.*DX/288.
      WT(2)=375.*DX/288.
      WT(3)=250.*DX/288.
      WT(4)=WT(3)
      WT(5)=WT(2)
C     MULTIPLY FA BY INTEGRATION WEIGHTS
      DO 130 M=1,NX
      M1=M
      M2=2*NX+2 - M
      MX=MOD(M-1,5) + 1
      FA(M1 ) = FA(M1 )*WT(MX)
      FA(M2 ) = FA(M2 )*WT(MX)
  130 CONTINUE
C     ADJUST WT FOR END POINT R+i0 IN INTEGRATION RULE
      FA(NX1) = FA(NX1)*0.5*WT(1)
C
      IF((2*LTR+1)*(2*NX+2).GT.IBUFF) THEN
        WRITE(6,9900)LTR
 9900   FORMAT('0NOT ENOUGH CIRCULAR BUFFER SPACE FOR DEUTERON WAVE,'
     1  ,' LTR=' , I3, ' IS TOO LARGE')
        STOP 'BUFFER ERROR IN UNBIND'
      ENDIF
C
      INC=0
      IY1=0
      FS1=-FS(1)
      DO 1020 I=1,JR
      IY2=NS(1)*LPLUS
      FS2=-FS(2)
      DO 1010 J=1,JS
      LPX=0
      LPT=1
      DO 1000 LL=1,LPLUS
      LK=LL
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      FJ2=FLOAT(LL-1)+FS2
C     CALCULATE PROTON WAVE extension into complex plane
      LM=LL+LPLUS
      IZ2=IY2+LK
      IFLG=1
C     -----------------------------------------------------
      CALL UNBIND(LL-1,NX,NI,FP(1),D(IZ2),H(LM)
     1,           ETA(2)*fk(2),cfk2,DX2,R2,IFLG)
C     -----------------------------------------------------
C
      M2=2*NX+2
      LP=1
      IF(LL.EQ.1) LP=MIN0(LTR+1,LPLUS)
      IF(LPX.GE.LPLUS) GO TO 160
      DO 150 L=1,LP
C     CALCULATE DEUTERON WAVE AND PLACE IN CIRCULAR BUFFER
      LM=LPX+1
      IZ1=IY1+LM
      IFLG=0
C     -----------------------------------------------------
      CALL UNBIND(LPX,NX,NI,BUFFR(LPT),D(IZ1),H(LM)
     1,           ETA(1)*fk(1),cfk1,DX1,R1,IFLG)
C     -----------------------------------------------------
      LPX=LPX+1
      LPTBL(LPX)=LPT-1
      LPT=LPT+(NX+1)
      IF(LPT+(NX+1).GT.IBUFF/2) LPT=1
  150 CONTINUE
  160 CONTINUE
C
C     MULTIPLY FP BY FA and place in FP
C
      LM=LL+LPLUS
      CP = cmplx(COS(S(LM)-S(LPLUS+1)), SIN(S(LM)-S(LPLUS+1)) )
      NX2=2*NX+1
      DO 300 M=1,NX2
      FP(M) = CP*FP(M)*FA(M)
  300 CONTINUE
c
      LPT1=LPTBL(LP)+1
      LPT2=LPT1+(NX+1)-1
C
      DO 400 LP=LP1,LP2,2
      FJ1=FLOAT(LP-1)+FS1
      IND=LK+INC
      SUM1=0.0
      SUM3=0.0
C     CONSTRUCT OUTGOING WAVE AMPLITUDE FROM D = ( EXP(2i*delta) - 1)/2i
      IZ1=IY1+LP
      Smat = cmplx(0.0, 2.0)*D(IZ1) + 1.0
      CD = cmplx(COS(S(LP)-S(1)), SIN(S(LP)-S(1)) )
      IF(FJ1.LT.0.0.OR.FJ2.LT.0.0) GO TO 331
C     NUMERICAL INTEGRATION FOR RADIAL MATRIX ELEMENTS
C     INTEGRATION IS +i*Dy FOR +y AND -i*Dy FOR -y
C     WITH LIMITS ON y OF 0 --> infinity HENCE THE + SIGN
C     WHEN COMBINING OUTGOING AND INGOING DEUTERON WAVES
C
      FE=1.0
      NX1=NX+1
      DO 330 M=1,NX1
      M1=NX+2 - M
      M2=NX   + M
      IX1=LPTBL(LP)+M1
c
c     Let k = k_d - k_p - k_n,   EZ = exp(-k*dy),  FE = exp(-k*y)
c     use  UD(-, R-iy) = conjg[ UD(+, R+iy) ]
c
C     MULTIPLY BY DEUTERON OUTGOING WAVE PART, UD(+)*FP(R+iy)*exp(-k*y)
      SUM1 = SUM1 + FP(M1)*      BUFFR(IX1) *FE
C     MULTIPLY BY DEUTERON INGOING  WAVE PART, UD(-)*FP(R-iy)*exp(-k*y)
      SUM3 = SUM3 + FP(M2)*conjg(BUFFR(IX1))*FE
      FE=FE*EZ
  330 CONTINUE
  331 CONTINUE
C     UD = (S*UD(+) + UD(-))/2  (divide by 2i  NOT! since  dz = idy)
      FLL(IND) = FLL(IND)+(Smat*SUM1+SUM3)*(0.5*CD*ANORM)
      LK=LK+LPLUS
  400 CONTINUE
 1000 CONTINUE
      INC=INC+INCR
      IY2=IY2+LPLUS
      FS1=FS1+1.0
 1010 CONTINUE
      IY1=IY1+LPLUS
      FS2=FS2+1.0
 1020 CONTINUE
C
c     End of unbound state section
c
      ENDIF
C
 1100 CONTINUE
C
C
      REWIND 4
      JR=IBF(7)+1
      JS=IBF(8)+1
      IF(NLTR.NE.1) GO TO 1105
      IF(JR+JS.EQ.2) GO TO 2000
 1105 CONTINUE
C
C     WRITE RADIAL MATRIX ELEMENTS ON TAPE 2
C
      JSX=JS
      IF(IBF(6).EQ.0) JSX=1
      DO 1120 I=1,JR
      IY=I-1
      IF(IBF(5).EQ.0) IY=0
      DO 1110 J=1,JS
      IZ=J-1
      IF(IBF(6).EQ.0) IZ=0
      INC=INCR*(JSX*IY+IZ)+1
      INDEX=INC-1+INCR
      WRITE(2)(FLL(II),II=INC,INDEX)
 1110 CONTINUE
 1120 CONTINUE
 2000 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE UNBIND(L,NX,NI,Y,D,H,ETAK,cfk2,DX,R,IFLG)
C
c     Subroutine to extend wave functions into complex plane for unbound
c     stripping cases.
c
c     L      Orbital angular momentum
c     nx     Number of points in integration
c     ni     Number of intermediate points
c     y      Storage for resulting wave function
c     d      Scattering amplitude = (exp(2i*delta) - 1)/(2i)
c     h      Hankel functions at Rmax
c     etak   Coulomb parameter * wave number
c     cfk2   Square of wave number
c     dx     Integration step size
c     r      Rmax on real axis
c     iflg   =0 for single side of real axis,  =1 for both sides
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      double complex Y(*), D, H, cfk2, Y1, Z1, Z0, Z
     1             , T1, T3, Smat
c
      DATA Yinit/1.0E+00/
C
      R2=R*R
      FL1=L*(L+1)
      DX12=DX*DX/12.
      EZ = EXP(-DX*SQRT(real(cfk2)))  
c     set starting points
      X = DX*FLOAT(NX*NI+NI+1)
      Z = cmplx(R, X)
      ETAK2=2.0*ETAK
      IF(real(cfk2) .EQ. 0.0) THEN
        Y0 = 1.0 - cmplx(0.0, float(L)) *DX/Z
      ELSE
        Y0 = Yinit
      ENDIF
      Z0 = Y0 * (1.0 + DX12*(FL1/Z**2 + ETAK2/Z - cfk2))
C
      Z = Z - cmplx(0.0D+0, DX)
      Y1 = Yinit
      Z1 = Y1 * (1.0 + DX12*(FL1/Z**2 + ETAK2/Z - cfk2))
 
C     INTEGRATE INWARDS FROM R + iy TO R - iy FOR IFLG = 1
C                       FROM R + iy TO R      FOR IFLG = 0
C     REMOVE EXP( k*y ) FACTOR SO THERE IS NO OVERFLOW
C
      if(iflg .eq. 0) then
        n2 = nx+1
      else
        n2 = 2*nx+1
      endif
c
      DO 200 I=1,N2
      DO 190 J=1,NI
      Z = Z - cmplx(0.0D+0, DX)
      T1=12.0*Y1
      T1=(T1 - 10.0*Z1 - Z0*EZ)*EZ
      Z0 = Z1
      Z1 = T1
      Y1 = Z1 / (1.0 + DX12*(FL1/Z**2 + ETAK2/Z - cfk2))
  190 CONTINUE
      Y(I)=Y1
  200 CONTINUE
C
      IF(IFLG .EQ. 0) then
C
C     NORMALIZE 'DEUTERON' WAVE TO OUTGOING HANKEL FUNCTION
C
C     Z0 is normalization constant to an
C     outgoing Hankel fct at Z = RMAX + i0
C     Z0 = H(+)/Y( R + i0)
      Z0 = H/Y(NX+1)
      DO 220 I=1,NX+1
      Y(I ) = Y(I )*Z0
  220 CONTINUE
c
      RETURN
      else
c
c     here for proton or neutron function
c
c     Construct an outgoing wave amplitude from H
C     H(1) IS G
C     H(2) IS F
C     Z0 is normalization constant to an
C     outgoing Hankel fct at Z = RMAX + i0
C     Z0   =  H(+)/Y(R + i0)
c
      Z0 = H/Y(NX+1)
c
C     D is scattering amplitude = (exp(2i*delta) -1)/(2i)
c
      Smat = cmplx(0.0, 2.0)*D + 1.0
      F1=1.0
      F2=EZ**(2*NI)
      NX1 = NX+1
      DO 250 I=1,NX1
      IK = NX+2 - I
      IL = NX   + I
C     T1 = exp(-k*X)*H(R-iy)
C     T3 = exp( k*X)*H(R+iy)
      T1 = Y(IL)*Z0
      T3 = Y(IK)*Z0*F1
c
C     u(R+iy) = exp(-k*X)*(S*H(R+iy) - H(R-iy))/(2*i)
      Y(IK) = (Smat*T3 - dconjg(T1) )*cmplx(0.0, -0.5)
C     u(R-iy) = exp(-k*X)*(S*H(R-iy) - H(R+iy))/(2i)
      Y(IL) = (Smat*T1 - dconjg(T3) )*cmplx(0.0, -0.5)
      F1=F1*F2
  250 CONTINUE
c
      endif
c
      RETURN
c
      END

c***********************************************************************
      SUBROUTINE BETAFN(FLL,D,LTR,JX,flfact,i_sym)
C
c     Subroutine to form scattering amplitude from radial matrix elements
c
c     FLL  Storage for radial matrix elements
c     D    Storage for scattering amplitudes
c     LTR  L  -transfer
c     JX   2*J-transfer
c     flfact normalizing factor
c     i_sym symmetry flags
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      double complex FLL(*), D(*)
      logical i_sym(2), igam1,igam2
C
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
c
      ISB=IS(3)
      JR=NS(1)
      JS=NS(2)
c
c     Set up spin for photon
c
      if(fm(1).eq.0.0.and.is(1).eq.2) then
        isi=0
        igam1=.true.
      else
        isi=is(1)
        igam1=.false.
      endif
      if(fm(2).eq.0.0.and.is(2).eq.2) then
        isf=0
        igam2=.true.
      else
        isf=is(2)
        igam2=.false.
      endif
c
      MPLUS=JX/2+1
      IFACT=MPLUS*JR*JS
      J2K=(1.0+PHASEF(NS(2)))/2.0
      M2K=JX-MPLUS-MPLUS+2
      LX=LTR+LTR
      TEMP2=SQRT(FLOAT((JX+1)*(IS(1)+1)))*flfact
      IF(FN.EQ.0.0) THEN
c       clear amplitude storage unless for coherent sum
        IND=2*LPLUS*IFACT
        DO 10 M=1,IND
        D(M)=0.0
   10   CONTINUE
      ENDIF
      IS1=-IS(1)
      DO 95 I=1,JR
      IS2=-IS(2)
      DO 90 J=1,JS
      IF(NLTR.NE.1) GO TO 14
      IF(JR*JS.EQ.1) GO TO 15
      IF(IBF(5)+IBF(6).EQ.0) GO TO 15
   14 CONTINUE
C
C     READ RADIAL MATRIX ELEMENTS FROM TAPE 2
C
      INCR = LPLUS*(LTR+1)
      READ (2)(FLL(INDEX),INDEX=1,INCR)
   15 continue
c     final L loop
      DO 80 LL=1,LPLUS
      lf=LL-1
      LLX=lf+lf
      JLX=LLX+IS2
      IF(JLX.LT.0) GO TO 75
      if(i_sym(2)) then
        if(phasef(lf).gt.0.0) then
          temp4=2.0*temp2
        else
          temp4=0.0
        endif
      else
        temp4=temp2
      endif      
      TEMP4=temp4*SQRT(FLOAT(JLX+1))*float(llx+1)
      if(igam2) then
         temp4=temp4*sqrt(float(ll)/(float(lf)+1.e-6))
     1              *(z(1)+za(1)*(-fm(1)/fma(1))**lf)
      endif
      LSTOR=lf*IFACT
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      LK=0
c     initial L loop
      DO 60 LP=LP1,LP2,2
      li=lp-1
      LPX=LP+LP-2
      JPX=LPX+IS1
      IF(JPX.GE.0) then
      if(i_sym(1)) then
        if(phasef(li).gt.0.0) then
          temp3=2.0*temp4
        else
          temp3=0.0
        endif
      else
        temp3=temp4
      endif      
      if(igam1) then
         temp3=temp3*sqrt(float(lp)/(float(li)+1.e-6))
     1              *(z(2)+za(1)*(-fm(2)/fma(2))**li)
      endif
      TEMP=temp3*SQRT(FLOAT(LPX+1))*PHASEF((LP-LL-LTR)/2)
     1    *VCC(LLX,LX,LPX,0,0)*WINEJ(LLX,isf,JLX,LX,ISB,JX,LPX,isi,JPX)
      INDEX=LK+LL
      KT=0
c      Initial state spins
      MSP=-IS(1)
      DO 57 MS1=1,JR
c      Final   state spins
      MS =-IS(2)
      DO 55 MS2=1,JS
      VCP=VCC(LPX,IS(1),JPX,0,MSP)
c
      DO 50 M=1,MPLUS
      MK=M+M-1
      MX=MK-1+M2K
      ML2=MSP-MX-MS
      ML=IABS(ML2/2)
      IF(ML.LE.lf) then
        IND=LSTOR+KT+M
        FACT=VCP*VCC(JLX,JX,JPX,MSP-MX,MX)*VCC(LLX,IS(2),JLX,ML2,MS)
     1      *SQRT(YXFCT(lf+ML,lf-ML))*TEMP
        D(IND)=D(IND)+FLL(INDEX)*FACT
      endif
   50 CONTINUE
      KT = KT+MPLUS
      MS =MS +2
   55 CONTINUE
      MSP=MSP+2
   57 CONTINUE
      endif
      LK=LK+LPLUS
   60 CONTINUE
   75 CONTINUE
   80 CONTINUE
      IS2=IS2+2
   90 CONTINUE
      IS1=IS1+2
   95 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE INSIG(D,PLM,JTR,FACTR)
c
c     Calculate inelastic cross sections and spin observables.
c
c     D		Final state amplitudes
c     PLM	Legendre polynomial storage
c     JTR	2*J_transfer
c     Factr	Normalization factor from BDWUCK4
c***********************************************************************
c
      parameter( ispc0 = 4010, ispc1 = 8000, ispc2 = 8000)
      implicit real*8(a-h,o-z)
      parameter(pi = 3.14159265, NX = 200, npol = 10)
      double complex D(*), SUM1, SUM(1000)
      logical i_open20, i_out20, i_20flag
c
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array2/Space2(ispc2)
c
      DIMENSION PLM(1000), Polz(npol)
      DIMENSION Sigplt(NX), Asyplt(NX)
      EQUIVALENCE (Space2(1), SUM(1)), (Space2(3601), Sigplt(1))
     1          , (Space2(3601+NX), Asyplt(1))
     2          , (ANGLE(1), THETAN),(ANGLE(2), THETA1)
     3          , (ANGLE(3), DTHETA), (SIG, Polz(1))
      data i_20flag / .true. /
C
      if(icon(17) .eq. 2 .or. icon(17) .eq. 3) then
        i_out20  = .true.
        if(i_20flag) then
          i_open20 = .true.
        endif  
        i_20flag = .false.
      else  
        i_open20 = .false.
        i_out20  = .false.
      endif
      NTHETA=THETAN
      NTHETA = min0(NTHETA, NX)
      asymax=0.0
      JR=NS(1)
      JS=NS(2)
      if(is(1).eq.2.and.fm(1).eq.0.0) then
c     initial state average factor for Gamma ray
        FACTR=sqrt(FACTR*3.0/2.0)
      else
        FACTA=sqrt(FACTR)
      endif
c
      M2K=(1.0-PHASEF(IS(3)))/2.0
      NPLUS=(JTR+IS(1)+IS(2))/2+1
      MPLUS=JTR/2+1
      IFACT = MPLUS*JR*JS
      WRITE(6,9000)
      IF(NTHETA.EQ.0) GO TO 230
      TotSig=0.0
c
      THETA=THETA1
      DO 120 NTH=1,NTHETA
      CALL LGNDR(PLM,NPLUS,LPLUS,THETA)
c
      Index1 = JS*JR*((JTR+1)/2)
      Index2 = Index1 + 1
      DO 100 M=1,MPLUS
      M2 = 2*(M-1)+M2K
      KT=M
      IS1=-IS(1)
c     Loop for initial spin substates
      DO 70 I=1,JR
      IS2=-IS(2)
c     Loop for final   spin substates
      DO 60 J=1,JS
      ML=-(M2+IS2-IS1)/2
      PHAS1 = 1.0
      IF(ML .lt. 0) PHAS1 = PHASEF(ML)
      PHAS2 = 1.0
      IF(ML .gt. 0) PHAS2 = PHASEF(ML)
      ML1 = IABS(ML)*LPLUS
      SUM1=0.0
c
      IND=KT
      DO 40 LL=1,LPLUS
      ML1 =ML1 +1
      SUM1 = SUM1+D(IND)*PLM(ML1)
C
C     CALCULATE TOTAL INELASTIC SIGMA
C
      IF(NTH .eq. 1) THEN
        L=LL-1
        ML = iabs(ML)
        if(ML.le.L) then
          FACT = conjg(D(IND))*D(IND)*YXFCT(L-ML,L+ML)/FLOAT(2*L+1)
          IF(M2 .NE. 0) FACT=FACT*2.0
          TotSig=TotSig+FACT
        endif
      ENDIF
      IND=IND+IFACT
   40 CONTINUE
        index1 = index1+1
        SUM(index1) = SUM1*PHAS1 *FACTA
      if(M2 .ne. 0) then
        index2 = index2-1
        SUM(index2) = SUM1*PHAS2 *FACTA
      endif
c	if(nth.eq.2) write(*,'(a,4i3, 1p4e12.4)')
c     1  ' Is2,Is1 M, ML :',is2,is1,M,ML,SUM(Index1),SUM(Index2)
      KT = KT+MPLUS
      IS2=IS2+2
   60 CONTINUE
      IS1=IS1+2
   70 CONTINUE
c
  100 CONTINUE
c
      Maxi = JTR + 1
      Max1 = 1
      CALL POLFCT(Max1,Maxi,JR,JS,Theta,Polz,SUM
     1  ,i_open20,i_out20,nth,ntheta,ALPHA,IDAT)
c
      WRITE(6,9001)THETA,(Polz(I), I=1,10),THETA
      SIGPLT(NTH)=SIG
      ASYPLT(NTH)=Polz(3)
      asymax=max(asymax,abs(asyplt(nth)))
      THETA=THETA+DTHETA
  120 CONTINUE
c
      TotSig =TotSig *4.0*pi*FACTA**2/float(JR)
      WRITE(6,9002)TotSig
c
      IF(ICON(13).ne.0) then
        WRITE(7,9011)ALPHA
        WRITE(7,9010)(SIGPLT(N),N=1,NTHETA)
      ENDIF
      IF(ICON(9).NE.0) then
        NTH=NTHETA
        WRITE(6,9011)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,ICON(9),SIGPLT,DTHETA,THETA1)
        NTEMP=-10
          if(IS(1).NE.0.and.asymax.gt.0.05) then
            WRITE(6,9013)ALPHA,(IDAT(I),I=1,3)
            CALL DWPLOT(NTH,NTEMP  ,ASYPLT,DTHETA,THETA1)
          endif
        ENDIF
  230 continue
      RETURN
c
 9000 FORMAT('0 Theta',' Inelsig,fm**2'
     1,   '      Polz','       Asy','       Ayy'  
     2,   '       A22','       A21','       A20'
     3,   '       T22','       T21','       T20'
     4,   '  Theta')
 9001 FORMAT(0PF7.2,1PE14.4,0P9F10.4,0PF7.2)
 9002 FORMAT('0Tot-sig',1PE13.4)
 9010 FORMAT(1PE10.4,5(4X,1PE10.4))
 9011 FORMAT('1Inelastic  ',15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9013 FORMAT('1Asymmetry  ',15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

