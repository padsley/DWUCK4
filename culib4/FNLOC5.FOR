
      SUBROUTINE FNLOC5(U,V,W,PNLOC,FK2,FK,ETA,RC,DR,KT)
c
c     W       Central optical potential
c     V       Spin orbit potential, for Darwin form
c     W       Function for correction to be applied to
c     PNLOC   Non-locality parameter
c             Positive value gives Gaussian form
c             Negative value gives Darwin form
c     FK2     Wave numbre squared
c     FK      Wave number
c     ETA     Coulomb parameter
c     RC      Coulomb radius
c     DR      Radial increment
c     KT      Number of points
c
      DIMENSION U(800),V(800),W(800)
C
C     NON LOCALITY CORRECTION FACTOR FOR DISTORTED WAVES
C
      FACP=PNLOC**2/8.0
      IF(FACP.EQ.0.0) GO TO 410
C
      WRITE(6,9900)PNLOC
      IF(PNLOC.GT.0.0) THEN
      R=0.0
      ETAK=2.0*ETA*FK
      ELSE
      R1=DR*FLOAT(KT+1)
      ENDIF
      CTEMP1=0.0
      CTEMP2=0.0
      T1=0.0
      T2=0.0
C
      DO 400 M=1,KT
      IF(PNLOC.GT.0.0) THEN
C     PNLOC POSITIVE GIVES USUAL NON-LOCAL FACTOR WITH GAUSSIAN SHAPE
      MK=M+M-1
      R=R+DR
      IF(R.GT.RC) THEN
      VCOUL=ETAK/R
      ELSE
      VCOUL=0.5*ETAK*(3.0-(R/RC)**2)/RC
      ENDIF
      CTEMP1=-FACP*U(MK  )+FACP*(FK2-VCOUL)
      CTEMP2=-FACP*U(MK+1)
      T1=U(MK  )
      T2=U(MK+1)
      ELSE
C     PNLOC NEGATIVE GIVES DIRAC TYPE DARWIN FACTOR BASED ON SPIN ORBIT
      MK=KT+KT-M-M+1
      R2=R1
      R1=R1-DR
      CTEMP1=CTEMP1-(R2*T1+R1*V(MK  ))*(DR*0.25)
      CTEMP2=CTEMP2-(R2*T2+R1*V(MK+1))*(DR*0.25)
      T1=V(MK  )
      T2=V(MK+1)
      ENDIF
C
      FACT=EXP(CTEMP1)
      UT1=FACT*COS(CTEMP2)
      UT2=FACT*SIN(CTEMP2)
      UF1   =UT1*W(MK  )-UT2*W(MK+1)
      UF2   =UT1*W(MK+1)+UT2*W(MK  )
      W(MK  )=UF1
      W(MK+1)=UF2
C	WRITE(20,7777)R1,T1,T2,UT1,UT2
C 7777	FORMAT(' R =',F6.2,'  POT =',1P2E12.4,'  FACTOR =',1P2E12.4)
  400 CONTINUE
  410 CONTINUE
      RETURN
 9900 FORMAT('0Non-locality correction is applied, parameter =',f8.3)
      END