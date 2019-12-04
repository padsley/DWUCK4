c$debug
c***********************************************************************
      SUBROUTINE ELSIG(dtemp,d,plm,SigPlt,angle,fk,eta,rsig,alpha
     1                ,idat,is,icon,lplus,isym)
c
c     Subroutine to print out elastic scattering cross sections
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 SUM(25),dtemp(*),d(*)
      logical isym(2), i_open20, i_out20
      DIMENSION plm(1000),SigPlt(200)
     1,angle(3),fk(2),eta(2),rsig(2),alpha(15),idat(6),is(2),icon(20)
     2,POL(12)

      DATA NX, N2, Pi/200, 1, 3.141592654/
c
      i_open20 = .false.
      i_out20  = .false.
C
      NTHETA=angle(1)
      theta1=angle(2)
      dtheta=angle(3)
      IF(NTHETA.EQ.0) return
c
      ICNT=0
      IC1=0
      DO 200 N=1,N2
      ISS=IS(N)+1
      IS2=IS(N)
      IQ=LPLUS*ISS**2
C     Clear Amplitude Storage
      DO 3 IP=1,IQ
      DTEMP(IP)=0.0
    3 CONTINUE
C
C     CONSTRUCT TABLE OF SCATTERING AMPLITUDES WITH ANG-MOM COEFFICENTS
C
      DO 25 LL=1,LPLUS
      L=LL-1
      LK=LL+LL-1
      L2=LK-1
      FL=LK
      J1=L2-IS2
      J2=L2+IS2
      IQ=IC1+LL
c
      DO 20 JJ2=J1,J2,2
      IP=LL
      DO 15 MS2=-IS2,IS2,2
      VCI=VCC(L2,IS2,JJ2,0,MS2)*FL
      DO 14 MS1=-IS2,IS2,2
      ML2=MS2-MS1
      ML=IABS(ML2/2)
      IF(ML .LE. L) THEN
	VCF=VCC(L2,IS2,JJ2,ML2,MS1)
        DTEMP(IP)=DTEMP(IP)+D(IQ)
     1     *VCI*VCF*SQRT(YXFCT(L+ML,L-ML))
      ENDIF
      IP=IP+LPLUS
   14 CONTINUE
   15 CONTINUE
       IQ=IQ+LPLUS
   20 CONTINUE
   25 CONTINUE
      IC1=IC1+LPLUS*ISS
C
C     CALCULATE ELASTIC CROSS SECTION
C
      THETA=THETA1
      if(ETA(N) .ne. 0.0) then
        WRITE(6,9027)
      else
        WRITE(6,9028)
      endif
      DO 100 NTH=1,NTHETA
      CALL LGNDR(PLM,ISS,LPLUS,THETA)
      Nss = ISS**2
      DO 30 I=1,Nss
      SUM(I)=0.0
   30 CONTINUE
      DO 40 LL=1,LPLUS
      IP=LL
      DO 35 I=1,ISS
      DO 34 J=1,ISS
      MP = iabs(J-I)*LPLUS+LL
      index = (I-1)*ISS + J
      Phas  = 1.0
      if(I .lt. J) Phas = phasef(I-J)
      SUM(index) = SUM(index) + DTEMP(IP)*PLM(MP)*Phas
      IP=IP+LPLUS
   34 CONTINUE
   35 CONTINUE
   40 CONTINUE
C
C     CALCULATE COULOMB AMPLITUDE
C
      if(theta.lt.0.0625
     1        .or. (isym(N) .and. abs(theta-180.).lt. 0.0625)) then
        CL=0.0
        SL=0.0
      else
        ARG  =THETA*(Pi/360.)
        S    =SIN(ARG)**2
        ARG  =ETA(N)*LOG(S)
        FFACT=ETA(N)/(2.0*FK(N)*S)
        CL   =   COS(ARG)*FFACT
        SL   =   SIN(ARG)*FFACT
      if(isym(N) .and. is(N).eq.0) then
        ARG  =THETA*(Pi/360.)
        C    =COS(ARG)**2
        ARG  =ETA(N)*LOG(C)
        GFACT=ETA(N)/(2.0*FK(N)*C)
        CL   =CL+COS(ARG)*GFACT
        SL   =SL+SIN(ARG)*GFACT
      endif
      endif
c     Add Coulomb amplitude
      DO 60 I=1,ISS
      index = (I-1)*ISS + I
        SUM(index) = SUM(index) - cmplx(CL, -SL)
   60 CONTINUE
c
      CALL POLFCT(1,1,ISS,ISS, theta,POL,SUM
     1          ,i_open20,i_out20, nth,nthetra, ALPHA,IDAT)
      if(eta(n).eq.0.0) then
        Ratio = 1.0
        SigPlt(NTH) = Pol(1)
      else
        if(theta.lt.0.0625
     1          .or. (isym(N) .and. abs(theta-180.).lt.0.0625)) then
          Ratio  = 1.0
          Pol(1) = 0.0
        else
          Ratio = Pol(1)/(CL**2+SL**2)
        endif
        SigPlt(nth) = Ratio
      endif
      cs = cos(theta*(Pi/180.))
      if(ETA(N) .ne. 0.0) then
        WRITE(6,9030)THETA,Ratio,(POL(I),I=1,7),THETA,cs
      else
        WRITE(6,9031)THETA,      (POL(I),I=1,7),THETA,cs
      endif
   90 CONTINUE
      THETA=THETA+DTHETA
  100 CONTINUE
      WRITE(6,9002)RSIG(N)
      NTH=MIN0(NTHETA,NX)
c
      IF(ICON(6).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,ICON(6),SigPlt,DTHETA,THETA1)
      endif
  200 continue
      RETURN
c
 9002 FORMAT('0REACSIG  ',1PE13.4)
 9027 FORMAT('0  Theta    Sig(1)/Coul     Sigma(1)  ','    Pol   '
     1       ,'    Asy   ','    Ayy   ','    A22   ','    A21   '
     2       ,'    A20   ','   Theta  ','Cos(Theta)')
 9028 FORMAT('0  Theta     Sigma(1)  ','    Pol   '
     1       ,'    Asy   ','    Ayy   ','    A22   ','    A21   '
     2       ,'    A20   ','   Theta  ','Cos(Theta)')
 9030 FORMAT(0PF8.2,2X,1P2E13.4,0P6F10.4,0PF10.2,0pf10.3)
 9031 FORMAT(0PF8.2,2X,1P1E13.4,0P6F10.4,0PF10.2,0pf10.3)
 9999 FORMAT('1',15A4,I4,2('/',I2.2),I4,2('.',I2.2))
      END