
      SUBROUTINE SLATR (KT,KMAX,DRF,VB,MINL,FMU,ICODE)
c
c     Calculates expansion of potentials for Cathen
c
      DIMENSION VB(800)
C
C     YUKAWA SLATER EXPANSION
C
      R=0.0
      IF(ICODE.EQ.2) GO TO 101
C
C     HANKEL FUNCTION*EXP(+FMU*R)
C
      F1=1.0
      F2=EXP(-FMU*DRF)
      DO 100 M=1,KMAX
      R=R+DRF
      X=FMU*R
      F1=F1*F2
      AZ=(1.0-F1*F1)/(2.0*X)
      B2=1.0/X
      B1=B2*(1.0+B2)
      FL=-1.0
      DO 50 LL=1,MINL
      B3=FL*B2/X+B1
      B1=B2
      B2=B3
      FL=FL+2.0
   50 CONTINUE
      VB(M+KMAX)=B2
C
C     BESSEL FUNCTION*EXP(-FMU*R)
C
      MAX=2.0*X+10.0
      MAX=MAX0(MAX,MINL+5)
      A3=0.0
      A2=1.0
      FL=MAX+MAX+3
      DO 70 LL=1,MAX
      N=MAX+1-LL
      FL=FL-2.0
      A1=A3+FL*A2/X
      IF(N.LT.MINL) GO TO 69
      IF(ABS(A1).LT.1.0E+20) GO TO 69
      A1=A1*1.0E-20
      A2=A2*1.0E-20
   69 CONTINUE
      IF(N.EQ.MINL) TEMP=A1
      A3=A2
      A2=A1
   70 CONTINUE
      VB(M    )=TEMP*AZ/A1
  100 CONTINUE
      RETURN
C
C     COULOMB SLATER EXPANSION
C
  101 CONTINUE
      FL=MINL+MINL-1
      DO 200 M=1,KMAX
      R=R+DRF
      A2=1.0
      DO 105 LL=1,MINL
      A2=A2*R
  105 CONTINUE
      VB(M    )=A2/(R*FL)
      VB(M+KMAX)=1.0/A2
  200 CONTINUE
      RETURN
      END

      SUBROUTINE RADIN(KT,KMAX,DRF,FMU,VB,UB,UC,SL,OPT,SI,KMT,SK)
c
c     Calculates slater integrals for Cathen
c
      DIMENSION UB(400),UC(400),VB(800),SI(800),DG(2)
      DATA XRHO/80./
c      FLOAT(III)=DFLOAT(III)
C
      KM1=KMAX
      KM2=KMAX*2
      KM3=KMAX*3
      DG(1)=2.0*DRF/3.0
      DG(2)=2.0*DG(1)
      F0=EXP( FMU*DRF)
      MMAX=XRHO/FMU/DRF
      R2=0.0
      SUMA=0.0
      SUMB=0.0
      SUMC=0.0
      SUMD=0.0
      DO 200 M=1,KT
      MK=M+M-1
      R2=R2+DRF
      F2=EXP(-FLOAT(MIN0(MMAX,M)-1)*FMU*DRF)
      SLL=0.0
      MX=1
      R1=0.0
      DO 100 MM=1,KT
      R1=R1+DRF
      MX=3-MX
      if(iabs(m-mm).lt.mmax) then
        if(m.lt.mm) then
          F2=F2/F0
          TEMP=VB(M)*VB(MM+KM3)*F2
        elseif(m.gt.mm) then
          TEMP=VB(MM+KM2)*VB(M+KM1)*F2
          F2=F2*F0
        else
          TEMP=0.5*(VB(M)*VB(M+KM3)+VB(M+KM2)*VB(M+KM1))
        endif
        F1=DG(MX)*UB(MM)*UC(MM)*R1**2
        SLL=SLL+F1*TEMP
      endif
  100 CONTINUE
      SLL=SLL*OPT
      R22=R2**2
      IF(M.EQ.KMT) SK=SLL
      SI(MK)=SI(MK)+SLL
      SUMA=SUMA+UB(M)*UC(M)*R22
      SUMB=SUMB+UB(M)*UC(M)*R22**2
      SUMC=SUMC+SLL*R22
      SUMD=SUMD+SLL*R22**2
  200 CONTINUE
      SUMA=SUMA*DRF
      SUMB=SUMB*DRF
      SUMC=SUMC*DRF
      SUMD=SUMD*DRF
      SL=SUMC
      WRITE(6,9100)SUMA,SUMB,SUMC,SUMD
      RETURN
 9100 FORMAT(13H0        J0 =,F11.4,7H   J1 =,F11.4,7H   K0 =,F11.4
     1 ,7H   K1 =,F11.4)
      END