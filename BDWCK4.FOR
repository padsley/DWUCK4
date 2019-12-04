c***********************************************************************
      SUBROUTINE BDWUCK4
c
c     Control program for second part of DWUCK4 which integrates the 
c     distorted waves, calculates radial integrals, inelastic
c     amplitudes and the crossections.
c     
c***********************************************************************
c
      parameter(ispc0 = 4010, ispc1 = 8000, ispc2 = 8000)
      IMPLICIT REAL*8(A-H,O-Z)
      logical i_sym(2)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array0/space0(ispc0)
      Common/array1/space1(ispc1)
      Common/array2/space2(ispc2)

      DIMENSION D(8000),sigplt(200,2),plm(2000), dtemp(2000), F(800)
     1,         FLL(8000), UB(800), FF(10)
      Equivalence (space1(   1), D), (space2(  1), FLL)
     1,           (space0(   1), F), (space0(   1), plm(1))
     2,           (space0(1601), sigplt(1,1)), (space0(3201), UB)
     3,           (space0(4001), FF) ,(space1(1001), dtemp)
      DATA FOURPI/12.5663706/
c
      do 40 i=1,2
c
      if(fm(i).eq.fma(i).and.z(i).eq.za(i)) then
        i_sym(i)=.true.
      else
        i_sym(i)=.false.
      endif
   40 continue
C
      IBUF=IBF(1)
      write(*,'(a)')'   Subroutine INTEG  entered'
      CALL INTEG4(i_sym)
C
C     PRINT OUT ELASTIC CROSS-SECTIONS 
C
      IF(ICON(6).ne.0.or.ICON(16).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        WRITE(6,9904)
C
        CALL ELSIG(dtemp,d,plm,sigplt,angle,fk,eta,rsig,alpha
     1            ,idat,is,icon,lplus,i_sym)
C
      else
        WRITE(6,9002)(RSIG(N),N=1,2)
      endif
      IF(ICON(12).NE.0) then
c
c     print out radial wave functions
c
c     skip over form factors
      DO 72 M=1,NLTR
      READ(4)
   72 CONTINUE
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
C
      CALL TAPED
C
      endif
c
      RZ=FLOAT(KZ)*DRF
      RMAX=FLOAT(K)*DRF
      JR=NS(1)
      JS=NS(2)
      DO 101 II=1,NLTR
      LTR=LTRT(II)
      MPLUS=LTR+1
C
C     SPACE TO PROPER FORM FACTOR FOR THIS L-TRANSFER
C
      DO 70 IJ=1,NLTR
      IF(IJ.EQ.II) then
        READ (4)UB,FF
      else
        READ (4)
      endif
   70 continue
C
      write(*,'(a)')'   Subroutine RADINT entered'
      CALL RADINT(F, FLL, UB, FF, LTR)
C
      IF(ICON(7).ne.0) then
C
C     WRITE( RADIAL MATRIX ELEMENTS IF ICON(7).NE.0
C
      INCR=LPL2*LTR
      INC=1
      IS1=-IS(1)

      DO 100 I=1,JR
      IS2=-IS(2)
      DO  98 J=1,JS
      WRITE(6,9999)ALPHA,(IDAT(M),M=1,3)
      WRITE(6,9900)IS1,IS2
      WRITE(6,9901)
      DO 90 LL=1,LPLUS
      LM=LL-1
      IND=INC+INCR
      WRITE(6,9902)LM,(FLL(IT),FLL(IT+1),IT=INC,IND,LPL2)
      INC=INC+2
   90 CONTINUE
      INC=INC+INCR
      IS2=IS2+2
   98 CONTINUE
      IS1=IS1+2
  100 CONTINUE
      endif
  101 CONTINUE
      END FILE 2
      REWIND 2
C
C     RESTORE SPIN STORAGE
C
      IS(1)=IBF(7)
      FS(1)=FLOAT(IBF(7))/2.0
      NS(1)=IBF(7)+1
      JR=NS(1)
      IS(2)=IBF(8)
      FS(2)=FLOAT(IBF(8))/2.0
      NS(2)=IBF(8)+1
      JS=NS(2)
      DO 300 II=1,NLTR
      LTR=LTRT(II)
      JTR=JTRT(II)
      IS(3)=ISTRT(II)
C
C     CALCULATE NORMALIZATION FACTOR
C
c     fact normalizes the cross section
c     flfact normalizes the amplitudes
c
      flfact=1.0
C
c     Photo-capture
      if    (fm(2).eq.0.0.and.is(2).eq.2) then
c       (p,gamma) reaction normalization
        fact=     2.0*fmu(1)*amu/(hbarc*fk(1))**2
     1            *chsq/fk(2)**2
        flfact=sqrt(float(2*ltr+1))
c
c     Photo-disintegration
      elseif(fm(1).eq.0.0.and.is(1).eq.2) then
c       (gamma,p) reaction normalization 
        fact=     2.0*fmu(2)*amu/(hbarc*fk(2))**2
     1            *chsq/fk(1)**2
      else
          FACT=     2.0*FMU(1)/(HBARC*FK(1))**2
     1             *2.0*FMU(2)/(HBARC*FK(2))**2
     2             *AMU**2/FOURPI
        if(abs(fm(1)-fm(2)).GT.0.1) then
C         Stripping normalization factors
          flfact=100.0*SQRT(FLOAT(2*LTR+1)/FLOAT(JTR+1))
          fact=fact*float(jtr+1)
        endif
      endif
c
      FACT=FACT*FK(2)/FK(1)
      FN=(II-1)*ICON(3)*(ICON(3)-1)
C
      CALL BETAFN(FLL, D, LTR,JTR,flfact,i_sym)
C
      MPLUS=JTR/2+1
      IF(icon(3).le.1.or.(icon(3).eq.2.and.II.eq.NLTR)) then
      IF(ICON(8).ne.0) then
C
C     WRITE BETA TABLES IF ICON(8).NE.0
C
      I1=MPLUS+MPLUS
      IFACT=I1     *NS(1)*NS(2)
      KT=1
      IS1=-IS(1)
      DO 120 I=1,JR
      IS2=-IS(2)
      DO 116 J=1,JS
      WRITE(6,9999)ALPHA,(IDAT(M),M=1,3)
      WRITE(6,9905)IS2,IS1
      WRITE(6,9903)
      DO 114 LL=1,LPLUS
      LM=LL-1
      IND=KT+LM*IFACT
      INDF=IND+I1-1
      WRITE(6,9902)LM,(D(INDEX),INDEX=IND,INDF)
  114 CONTINUE
      KT=KT+I1
      IS2=IS2+2
  116 CONTINUE
      IS1=IS1+2
  120 CONTINUE
      endif
C
C     SET HEADINGS FOR INELASTIC SIGMA
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9507)DRF,RZ,RMAX,VCE,FS(1)
      WRITE(6,9508)FNRNG,PNLOC,FS(2)
      TEMP=ECM(2)-ECM(1)
      do 280 jj=1,nltr
      if((icon(3).le.1.and.jj.eq.ii).or.
     1   (icon(3).eq.2.and.ii.eq.nltr)) then
        WRITE(6,9505)LTRT(jj),JTRT(jj),IS(3),TEMP
      endif
  280 continue
C
        write(*,'(a)')'   Subroutine INSIG  entered'
        CALL INSIG(D, plm, JTR,FACT)
C
      endif
  300 CONTINUE
      REWIND 2
      RETURN

 9002 FORMAT(1H0,'REACSIG 1',1PE14.5,45X,'REACSIG 2',1PE14.5)
 9505 FORMAT(18H ANG MOM TRANSFER ,9H   LTR  =,I4,14H.       2*JTR=,I4
     1,14H.       2*STR=,I4,14H.       Q    =,F9.4)
 9507 FORMAT(18H0BASIC DATA       ,9H   DR   =,F9.4,9H   RMIN =,F9.4
     1,                            9H   RMAX =,F9.4,9H   COUEX=,F9.4
     2,                            9H     IS1=,F9.4)
 9508 FORMAT(18X,9H   FNRNG=,F9.4,9H   PNLC1=,F9.4,9H   PNLC2=,F9.4
     1,                           9H   PNLC3=,F9.4,9H     IS2=,F9.4)
 9900 FORMAT(24H0 RADIAL MATRIX ELEMENTS ,9H,  J1=L1+,I2,2H/2
     1,                                   9H,  J2=L2+,I2,2H/2 )
 9901 FORMAT('0 L2  F(L2,/L2-LTR/  )    F(L2,/L2-LTR/+2)'
     1,        '    F(L2,/L2-LTR/+4)    F(L2,/L2-LTR/+6)')
 9902 FORMAT(I4,1P10E10.3/(4X,1P10E10.3))
 9903 FORMAT('0 L2  BETA(L2,0)          BETA(L2,1)          BETA(L2,2)'
     1,  '          BETA(L2,3)          BETA(L2,4)')
 9904 FORMAT(1H0,32HELASTIC SCATTERING CROSS-SECTION )
 9905 FORMAT(1H0,8H    MS2=,I2,2H/2,8H    MS1=,I2,2H/2)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

c***********************************************************************
      SUBROUTINE INTEG4(isym)
C
c     Subroutine to integrate radial differential equations for the
c     distorted waves.
c
c***********************************************************************
      parameter(ispc0 = 4010, ispc1 = 8000, ispc2 = 8000)
      IMPLICIT REAL*8(A-H,O-Z)
      logical isym(2)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array0/space0(ispc0)
      Common/array1/space1(ispc1)
      Common/array2/space2(ispc2)
      DIMENSION U(800,2), V(800,2)
     1         ,F (400),FP(400),G (400),GP(400),S(400)
     2         ,F1(800),F2(800),Q1(800),Q2(800),C(800)
     3         ,D (800),X(800),DTEMP(3200)
     4         ,E (12),Q(4),A(2),B(2),CTEMP(2),LM(6)
     5         ,DRR2(2),DR2(2),R(2)
      EQUIVALENCE (space0(   1), U), (space0(1601), V)
     1           ,(U,DTEMP)
     2           ,(space1(   1),Q1),(space1( 801),Q2),(space1(1601),F1)
     3           ,(space1(2401),F2),(space1( 801),X ),(space1(   1),D )
     4           ,(DTEMP(   1),F ),(DTEMP( 401),FP),(DTEMP( 801),G )
     5           ,(DTEMP(1201),GP),(DTEMP(1601),S ),(DTEMP(2001),C )
      DATA ETA3/10.E+00/
C
      IWORD=0
      JT=NS(1)+NS(2)
      NP=LPL2*JT
      I=0
      DO 30 N=1,2
      DR2(N)=DR(N)**2/12.0
      R(N)=0.0
      JS=NS(N)
      DO 29 ISS=1,JS
      I=I+1
      LM(I)=0
   29 CONTINUE
   30 CONTINUE
      DO 40 IQ=1,NP
      F1(IQ)=0.0
      F2(IQ)=0.0
      Q1(IQ)=0.0
      Q2(IQ)=0.0
   40 CONTINUE
C
      DO 100 M=1,K
      MK=M+M-1
      IX=0
      I=0
      DO 90 N=1,2
      R(N)=R(N)+DR(N)
      DRR2(N)=DR2(N)/R(N)**2
      Q(1)=1.0+DR2(N)*U(MK  ,N)
      Q(2)=    DR2(N)*U(MK+1,N)
      LTEMP=2.0*FK(N)*R(N)+ETA3
      LTEMP=MIN0(LTEMP,LPLUS)
      FI=-FS(N)
      JS=NS(N)
      SFACT=FS(N)**2+FS(N)
      DO 89 ISS=1,JS
      I=I+1
      FL=0.0
      DO 80 LL=1,LPLUS
      FJ=FL+FI
      IX1=IX+LL+LL-1
      FLFACT=FL**2+FL
      FACT=DR2(N)*(FJ**2+FJ-FLFACT-SFACT)*0.5
      Q(3     )=Q(1)+FACT*V(MK  ,N)-DRR2(N)*FLFACT
      Q(4     )=Q(2)+FACT*V(MK+1,N)
      IF(LL.LE.LM(I)) GO TO 70
      IF(LTEMP.LT.LL) GO TO 72
      LM(I)=LM(I)+1
      IF(FJ-ABS(FL-FS(N)).LT.0.0) GO TO 72
c       calculate approximate starting value
        f2(ix1  )=1.0
        do 50 ii=1,ll
        f2(ix1  )=f2(ix1  )*(fk(n)*r(n))/float(2*ii-1)
   50   continue
      F2(IX1+1)=0.0
      Q2(IX1  )=Q(3)*f2(ix1  )
      Q2(IX1+1)=Q(4)*f2(ix1  )
C
C     EVALUATE Q AT ORIGIN FOR L=1
C
      IF(LL.EQ.2) Q1(IX+3)=-f2(ix1  )/6.0
      GO TO 72
   70 CONTINUE
c
c     Step equations forward by dr(n) via Numerov-Fox-Goodwin-Milne method
c
      CTEMP(1)=12.*F2(IX1  )-10.*Q2(IX1  )-Q1(IX1  )
      CTEMP(2)=12.*F2(IX1+1)-10.*Q2(IX1+1)-Q1(IX1+1)
      F1(IX1  )=F2(IX1  )
      F1(IX1+1)=F2(IX1+1)
      DET=Q(3)**2+Q(4)**2
      F2(IX1  )=(CTEMP(1)*Q(3     )+CTEMP(2)*Q(4     ))/DET
      F2(IX1+1)=(CTEMP(2)*Q(3     )-CTEMP(1)*Q(4     ))/DET
      Q1(IX1  )=Q2(IX1  )
      Q1(IX1+1)=Q2(IX1+1)
      Q2(IX1  )=CTEMP(1)
      Q2(IX1+1)=CTEMP(2)
   72 CONTINUE
      FL=FL+1.0
   80 CONTINUE
      FI=FI+1.0
      IX=IX+LPL2
   89 CONTINUE
   90 CONTINUE
C
C     WRITE RADIAL WAVE FUNCTIONS ON TAPE 4
C
      WRITE(4)(F2(J),J=1,NP)
  100 CONTINUE
      LX=1
      drrc = 0.1
      DO 120 N=1,2
      R2=FK(N)*R(N)
      R1=R2-DR(N)*FK(N)
      CALL COU(R1,R2,ETA(N),LPLUS,drrc,F(LX),FP(LX),G(LX),GP(LX),S(LX))
      RSIG(N)=0.0
      LX=LX+LPLUS
  120 CONTINUE
C
      IF(ICON(5).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        WRITE(6,9600)
        WRITE(6,9601)
      endif
c
c     Match solutions to asymptotic form
c
      DO 300 LL=1,LPLUS
      FL=FLOAT(LL-1)
      LX=LL
      I=0
      IX1=LL+LL-1
      DO 200 N=1,2
      JS=NS(N)
      FI=-FS(N)
      ARG=S(LX)-S(LX-LL+1)
      Q(1)=COS(ARG)
      Q(2)=SIN(ARG)
      Q(3)=Q(1)**2-Q(2)**2
      Q(4)=2.0*Q(1)*Q(2)
      DO 199 ISS=1,JS
      FJ=FL+FI
      I=I+1
      DET=F(LX)*GP(LX)-FP(LX)*G(LX)
      A(1)=(F1(IX1  )*GP(LX)-F2(IX1  )*G (LX))/DET
      A(2)=(F1(IX1+1)*GP(LX)-F2(IX1+1)*G (LX))/DET
      B(1)=(F2(IX1  )*F (LX)-F1(IX1  )*FP(LX))/DET
      B(2)=(F2(IX1+1)*F (LX)-F1(IX1+1)*FP(LX))/DET
      IF(LL.LE.LM(I).and.FJ-ABS(FL-FS(N)).ge.0.0) then
        DET=(A(1)+B(2))**2+(A(2)-B(1))**2
        CTEMP(1)=(A(1)+B(2))/DET
        CTEMP(2)=(B(1)-A(2))/DET
      else
          CTEMP(1)=0.0
          CTEMP(2)=0.0
      endif
C
C     C=NORMALIZATION CONSTANTS
C
      C(IX1  )=Q(1)*CTEMP(1)-Q(2)*CTEMP(2)
      C(IX1+1)=Q(1)*CTEMP(2)+Q(2)*CTEMP(1)
C
C     E=PARTIAL WAVE SCATTERING AMPLITUDES
C
      E(2*I-1)=B(1)*CTEMP(1)-B(2)*CTEMP(2)
      E(2*I  )=B(1)*CTEMP(2)+B(2)*CTEMP(1)
      T1 = E(2*I-1)
      T2 = E(2*I  )
      if(isym(N) .and. is(N).eq.0 ) then
        if(Phasef(lx-1).lt.0.0) then
        T1 = 0.0
        T2 = 0.0
        E(2*I-1) = 0.0
        E(2*I  ) = 0.0
        else
        T1 = 2.0*T1
        T2 = 2.0*T2
        endif
      endif
C
C     D=PARTIAL WAVE SCATTERING AMPLITUDES* COULOMB PHASE / WAVE NUMBER
C
      D(IX1  )=(Q(3)*T1 - Q(4)*T2)/FK(N)
      D(IX1+1)=(Q(3)*T2 + Q(4)*T1)/FK(N)
      X(IX1  )=E(2*I-1)
      X(IX1+1)=E(2*I  )
C
C     CALCULATE REACTION SIGMA
C
      T1 = E(2*I  ) - E(2*I-1)**2-E(2*I  )**2 
      if(isym(N) .and. is(N).eq. 0) T1 = 4.0*T1 
      RSIG(N)=RSIG(N)+(2.0*FJ+1.0)*T1
      FI=FI+1.0
      IX1=IX1+LPL2
  199 CONTINUE
      LX=LX+LPLUS
  200 CONTINUE
      IF(ICON(5).ne.0) then
C
C     WRITE ELASTIC PARTIAL WAVE SCATTERING AMPLITUDES
C
        IX=LL-1
        I1=1
        I2=NS(1)*2
        WRITE(6,9602)IX,(E(INDEX  ),INDEX=I1,I2)
        I1=I2+1
        I2=I2+NS(2)*2
        WRITE(6,9603)IX,(E(INDEX  ),INDEX=I1,I2)
      endif
  300 CONTINUE
      DO 310 N=1,2
      RSIG(N)=RSIG(N)*12.566371 /((2.0*FS(N)+1.0)*FK(N)**2)
  310 CONTINUE
C
C     WRITE NORMALIZATION CONSTANTS ON TAPE 4
C
      WRITE(4)(C(I),I=1,NP)
C
C     WRITE END POINT FUNCTIONS FOR UNBOUND STATE STRIPPING
C
C     ILINE IS SET IN ADWUCK FORM FACTOR LOOP
      IF(ILINE.EQ.0) GO TO 400
      WRITE(4)(X(I),I=1,NP),(GP(I),FP(I),I=1,LPL2),(S(I),I=1,LPL2)
  400 CONTINUE
      END FILE 4
      REWIND 4
      RETURN
 9600 FORMAT(1H0,35HPARTIAL WAVE SCATTERING AMPLITUDES )
 9601 FORMAT(4H   L,20H   REAL D1   IMAG D1 ,20H   REAL D2   IMAG D2
     1,             20H   REAL D3   IMAG D3
     2,   4X,4H   L,20H   REAL D1   IMAG D1 ,20H   REAL D2   IMAG D2
     3,             20H   REAL D3   IMAG D3 )
 9602 FORMAT(1H ,    I3,6F10.6)
 9603 FORMAT(1H+,68X,I3,6F10.6)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

c***********************************************************************
      SUBROUTINE TAPED
C
c     Subroutine writes out radial wave functions,
c     spaced icon(12) points apart
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

      DIMENSION C(800),FR(800),R(2),CTEMP(2)
      EQUIVALENCE  (space0(2001), C(1)), (space0( 801), FR(1))
c
      IK=ICON(12)
      JT=NS(1)+NS(2)
      NP=LPL2*JT
      R(1)=0.0
      R(2)=0.0
      IKTEMP=0
      DO 50 M=1,K
      IKTEMP=IKTEMP+1
C
C     READ IN DISTORTED WAVE RADIAL FUNCTIONS
C
      READ (4)(FR(J),J=1,NP)
      IF(IKTEMP.NE.IK) GO TO 50
      IKTEMP=0
      IX=0
      DO 40 N=1,2
      R(N)=R(N)+DR(N)*FLOAT(IK)
      JX=NS(N)
      DO 39 J=1,JX
      DO 30 LL=1,LPLUS
      LK=LL+LL-1
      IX1=LK+IX
      IX2=IX1+1
C
C     NORMALIZE RADIAL FUNCTIONS
C
      CTEMP(1)=FR(IX1)*C(IX1)-FR(IX2)*C(IX2)
      CTEMP(2)=FR(IX1)*C(IX2)+FR(IX2)*C(IX1)
      FR(IX1)=CTEMP(1)
      FR(IX2)=CTEMP(2)
   30 CONTINUE
      IX=IX+LPL2
   39 CONTINUE
   40 CONTINUE
      WRITE(6,9001)R(1),R(2)
      WRITE(6,9601)
      DO 45 LL=1,LPLUS
      LM=LL-1
      LK1=LM+LL
      LK2=LPL2*(NS(1)-1)+LK1
      WRITE(6,9602)LM,(FR(LK),FR(LK+1),LK=LK1,LK2,LPL2)
      LK1=LK2+LPL2
      LK2=LPL2*(NS(2)-1)+LK1
      WRITE(6,9603)LM,(FR(LK),FR(LK+1),LK=LK1,LK2,LPL2)
   45 CONTINUE
   50 CONTINUE
      REWIND 4
      RETURN
 9001 FORMAT(1H0,3HR1=,F8.4,57X,3HR2=,F8.4)
 9601 FORMAT(4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2   IMAG D2 ,20H
     1 REAL D3   IMAG D3 ,4X,4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2
     2   IMAG D2 ,20H  REAL D3   IMAG D3 )
 9602 FORMAT(1H ,    I3,6F10.7)
 9603 FORMAT(1H+,68X,I3,6F10.7)
      END

