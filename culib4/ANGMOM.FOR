
c***********************************************************************
      FUNCTION VCC(JX1,JX2,JX3,MX1,MX2)
c
c     Clebsch-Gordan Coefficient Routine
c
c***********************************************************************
C      IMPLICIT REAL*8(A-H,O-Z)
c      EXTERNAL FACTOR
      COMMON/FACTRL/FACT(0:32)
C
      VCC=0.0
      J1=JX1
      J2=JX2
      J3=JX3
      M1=MX1
      M2=MX2
      IF(J1.LT.J2) GO TO 20
      IF(J3.LT.J2) GO TO 30
      ICNTR=0
      GO TO 40
   20 IF(J3.LT.J1) GO TO 30
      ICNTR=-1
      IT=J1
      J1=J2
      J2=IT
      IT=M1
      M1=M2
      M2=IT
      GO TO 40
   30 ICNTR=1
      IT=J2
      J2=J3
      J3=IT
      M2=-M1-M2
   40 CONTINUE
      JZ1=(J1+J2-J3)/2
      IF(JZ1.LT.0) GO TO 150
      JZ2=(J1+J3-J2)/2
      IF(JZ2.LT.0) GO TO 150
      JZ3=(J2+J3-J1)/2
      IF(JZ3.LT.0) GO TO 150
      IF(J1-IABS(M1).LT.0) GO TO 150
      IF(J2-IABS(M2).LT.0) GO TO 150
      IF(J3-IABS(M1+M2).LT.0) GO TO 150
      JT1=(J1-J3+M2)/2
      JT2=(J2-J3-M1)/2
      NUMIN=MAX0 (JT1,JT2,0)
      JT3=(J1-M1)/2
      JT4=(J2+M2)/2
      NUMAX=MIN0 (JT3,JT4,JZ1)
      JT5=(J2-M2)/2
      IF(NUMAX.LT.NUMIN) GO TO 150
      J4=J1/2
      J5=J3/2
      PHAS=PHASEF(NUMIN)
      DO 100 NU=NUMIN,NUMAX
      VCC=VCC+PHAS      *(YXFCT(JT3-NU,J4)*YXFCT(NU-JT2,J5))
     1/(FACT(JT4-NU)*FACT(NU-JT1)*FACT(JZ1-NU)*FACT(NU))
      PHAS=-PHAS
  100 CONTINUE
      FCTOR=YXFCT(J4,(J1+M1)/2)*YXFCT(J4,JT3)*YXFCT((J1+J2+J3)/2+1,JZ2)*
     1YXFCT(J5,(J3+M1+M2)/2)*YXFCT(J5,(J3-M1-M2)/2)*FACT(JZ1)*FACT(JZ3)*
     2FACT(JT4)*FACT(JT5)*FLOAT(J3+1)
      VCC=SQRT(FCTOR)*VCC
      IF(ICNTR)120,150,110
  110 VCC=VCC*SQRT(FLOAT(J2+1)/FLOAT(J3+1))*PHASEF(JT3)
      GO TO 150
  120 VCC=VCC*PHASEF(JZ1)
  150 RETURN
      END

c***********************************************************************
      FUNCTION PHASEF(N)
c
c     Calculates (-1)**N
c
c***********************************************************************
c      IMPLICIT REAL*8(A-H,O-Z)
c
      PHASEF=FLOAT(1-2*IABS(N-2*(N/2)))
      RETURN
      END

c***********************************************************************
      FUNCTION YXFCT(M,N)
c
C     COMPUTES N_FACTORIAL/M_FACTORIAL
c
c***********************************************************************
c      IMPLICIT REAL*8(A-H,O-Z)
c
      YXFCT=1.0
      NUMAX=M-N
      IF(NUMAX)30,100,20
   20 ICTRL=0
      FCTOR=N
      GO TO 40
   30 ICTRL=1
      NUMAX=-NUMAX
      FCTOR=M
   40 CONTINUE
      DO 50 NU=1,NUMAX
      FCTOR=FCTOR+1.0
      YXFCT=YXFCT*FCTOR
   50 CONTINUE
      IF(ICTRL.EQ.0) YXFCT=1.0/YXFCT
  100 RETURN
      END

c***********************************************************************
      FUNCTION RACAH(J1,J2,J3,J4,J5,J6)
c
c     Calculates Racah coefficients
c     Spins are input as twice the spin
c
c***********************************************************************
c
c      IMPLICIT REAL*8(A-H,O-Z)
      logical jy_big
c
c      EXTERNAL FACTOR
      COMMON/FACTRL/FACT(0:32)
c
      RACAH=0.0
      Z1=DELR(J1,J2,J5)
      IF(Z1.EQ.0.0) GO TO 90
      Z1=DELR(J3,J4,J5)*Z1
      IF(Z1.EQ.0.0) GO TO 90
      Z2=DELR(J1,J3,J6)
      IF(Z2.EQ.0.0) GO TO 90
      Z2=DELR(J2,J4,J6)*Z2
      IF(Z2.EQ.0.0) GO TO 90
      Z1=SQRT(Z1/Z2)*Z2
      JT1=(J1+J2+J5)/2
      JT2=(J3+J4+J5)/2
      JT3=(J1+J3+J6)/2
      JT4=(J2+J4+J6)/2
      JZ1=(J1+J2+J3+J4)/2
      JZ2=(J1+J4+J5+J6)/2
      JZ3=(J2+J3+J5+J6)/2
c
      NUMIN=MAX0 (JT1,JT2,JT3,JT4)
      NUMAX=MIN0 (JZ1,JZ2,JZ3)
      IF(NUMAX.ge.NUMIN) then
      if(NUMIN-JT4 .gt. JZ1-NUMIN) then
        jy_big = .true.
      else
        jy_big = .false.
      endif
      PHASE=PHASEF(NUMIN+JZ1)*Z1
      DO 80 NU=NUMIN,NUMAX
      JY1=NU-JT1
      JY2=NU-JT2
      JY3=NU-JT3
      JY4=NU-JT4
      JY5=JZ1-NU
      JY6=JZ2-NU
      JY7=JZ3-NU
      if(jy_big) then
        FACTR = YXFCT(JY4,NU+1)/FACT(JY5)
      else
        FACTR = YXFCT(JY5,NU+1)/FACT(JY4)
      endif
      RACAH=RACAH+PHASE*FACTR
     1  /(FACT(JY1)*FACT(JY2)*FACT(JY3)*FACT(JY6)*FACT(JY7))
      PHASE=-PHASE
   80 CONTINUE
      endif
   90 RETURN
      END

c***********************************************************************
      FUNCTION DELR(J1,J2,J3)
c
c     Triangular function
c     Used by Racah Function
c
c***********************************************************************
c      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FACTOR
      COMMON/FACTRL/FACT(0:32)
c
      JZ1=(J1+J2-J3)/2
      IF(JZ1.LT.0) GO TO 130
      JZ2=(J1-J2+J3)/2
      IF(JZ2.LT.0) GO TO 130
      JZ3=(J2+J3-J1)/2
      IF(JZ3.LT.0) GO TO 130
      JZ4=(J1+J2+J3)/2+1
      IF(JZ3.LT.JZ2) GO TO 80
      IF(JZ3.LT.JZ1) GO TO 70
      DELR=YXFCT(JZ4,JZ3)*FACT(JZ1)*FACT(JZ2)
      RETURN
   70 DELR=YXFCT(JZ4,JZ1)*FACT(JZ2)*FACT(JZ3)
      RETURN
   80 IF(JZ2.LT.JZ1) GO TO 70
      DELR=YXFCT(JZ4,JZ2)*FACT(JZ1)*FACT(JZ3)
      RETURN
  130 DELR=0.0
      RETURN
      END

c***********************************************************************
      FUNCTION WINEJ(J1,J2,J3,J4,J5,J6,J7,J8,J9)
c
c     NineJ Symbol Function
c
c***********************************************************************
c      IMPLICIT REAL*8(A-H,O-Z)
c
      WINEJ=0.0
      MUMIN=MAX0(IABS(J1-J9),IABS(J2-J6),IABS(J4-J8))
      MUMAX=MIN0(J1+J9,J2+J6,J4+J8)
      IF(MUMAX.LT.MUMIN) GO TO 40
       DO 20 MU=MUMIN,MUMAX,2
      PROD=RACAH(J1,J4,J9,J8,J7,MU)*RACAH(J2,J5,MU,J4,J8,J6)*
     1     RACAH(J9,MU,J3,J2,J1,J6)*FLOAT(MU+1)
      WINEJ=WINEJ+PROD
   20 CONTINUE
      WINEJ=WINEJ*PHASEF((J1+J3+J5+J8)/2+J2+J4+J9)
   40 RETURN
      END

      BLOCK DATA FACTOR
c
c     Factorial table
c
c      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FACTRL/FACT(0:32)
C
      DATA FACT/ 1.0000000000E+00, 1.0000000000E+00, 2.0000000000E+00
     1         , 6.0000000000E+00, 2.4000000000E+01, 1.2000000000E+02
     2         , 7.2000000000E+02, 5.0400000000E+03, 4.0320000000E+04
     3         , 3.6288000000E+05, 3.6288000000E+06, 3.9916800000E+07
     4         , 4.7900160000E+08, 6.2270208000E+09, 8.7178291200E+10
     5         , 1.3076743680E+12, 2.0922789888E+13, 3.5568742810E+14
     6         , 6.4023737057E+15, 1.2164510041E+17, 2.4329020082E+18
     7         , 5.1090942172E+19, 1.1240007278E+21, 2.5852016739E+22
     8         , 6.2044840173E+23, 1.5511210043E+25, 4.0329146113E+26
     9         , 1.0888869450E+28, 3.0488834461E+29, 8.8417619937E+30
     $         , 2.6525285981E+32, 8.2228386542E+33, 2.6313083693E+35/
C    $         , 8.6833176188D+36, 2.9523279904D+38, 1.0333147966D+40
C    $         , 3.7199332679D+41, 1.3763753091D+43, 5.2302261747D+44
C    $         , 2.0397882081D+46, 8.1591528325D+47, 3.3452526613D+49
C    $         , 1.4050061178D+51, 6.0415263063D+52, 2.6582715748D+54
C    $         , 1.1962222087D+56, 5.5026221598D+57, 2.5862324151D+59
C    $         , 1.2413915593D+61, 6.0828186403D+62, 3.0414093202D+64
C    $         , 1.5511187533D+66/
      END