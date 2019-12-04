      SUBROUTINE BIND(U,F,DR,RM,FNODE,FL,Kin,FK,ETA,V,E,FK2,ISW,IERR,D)
C
C     U       = POTENTIAL STORAGE
C        ODD  =             - SCALED PART OF V *(2M/HBAR**2)
C        EVEN = (ENERGY - NON SCALED PART OF V)*(2M/HBAR**2)
C     F       = OUTPUTTED WAVE FUNCTION
C     DR      = RADIAL INCREMENT
C     RM      = MATCHING RADIUS, INPUT IS DEFAULT, OUTPUT IS VALUE USED
C     FNODE   = NO. NODES IN WAVE FCT., EXCLUDING ORIGIN AND INFINITY
C     FL      = ORBITAL ANGULAR MOMENTUM
C     Kin     = MAX NO. OF POINTS + 2048*(MIN POINT)
C     FK      = |WAVE NUMBER| E.G. SIN(FK*R) OR EXP(-FK*R)
C     ETA     = COULOMB PARAMETER
C     V       = SCALE FACTOR FOR ODD POTENTIAL VALUES
C     E       = ENERGY
C     FK2     = E*(2*M/HBAR**2)
C     ISW     = SEARCH SWITCH
C        0    = SEARCH ON V
C        1    = SEARCH ON E
C        2    = NO SEARCH, FOR E>0 ONLY (PHASE SHIFT IS CALCULATED)
C        3    = Search on V, box b.c.
C        4    = Search on E. box b.c.
C     IERR    = ERROR FLAG, IS SET IF ERROR OCCURS
C     D       = LAST POINT DATA FOR UNBOUND STRIPPING
C     D(1),D(2) ARE REAL AND IMAG PARTS OF EXP(2I*DELTA)
C     D(3),D(4) ARE GL(RMAX) AND FL(RMAX)
C     D(5),D(6) are real and imaginary K**2
C     D(10)     is ISW
C
c      IMPLICIT REAL*8(A-H,O-Z)
c
      parameter (Pi = 3.141592654)
      logical dplus,dminus,turnpt
      CHARACTER*2 A
      DIMENSION U(800),F(400),T(8),D(10)
      common /ch3files/input,ioutp
      DATA ICMAX/32/
      DATA DLAMY,DETX,dinterp,dlammax/.001,.01,0.1,0.3/
      DATA XFACT, ELOG10, ONE/0.707, 2.302585, 1.0/
C
C     CHANGES TO ALLOW INTEGRATION OF HARD CORE POTENTIALS
C     THEY ARE KLUDGED TO ALLOW BACKWARD COMPATIBILITY
C     Kmin is initial point for hard core
C
	if(ioutp .eq. 0) then
	  ioutput = 6
	else
	  ioutput = ioutp
	endif
      Kmin = Kin/2048
      Kmax = Kin-Kmin*2048
      Kmin=Kmin+1
C
      CFACT=1.0
      dv=0.225*v
      v0=0.0
      v1=0.0
      d0 =-1.0e20
      d1 = 1.0e20
      dplus =.false.
      dminus=.false.
      IXNT=0
      IPRNT=0
      IF(IERR.LT.0) IPRNT=1
      IERR=0
      ICNT=ICMAX
      IF(ISW.EQ.2) ICNT=0
      FLP=FL*(FL+1.0)
      DR2=DR**2/12.0
      LL=FL+1.0
      ITER=0
C
C     CALCULATE OUTER BOUNDARY CONDITION
C
   10 CONTINUE
      RMAX=DR*FLOAT(Kmax)
      IF(FK2.GT.0.0) THEN
C     POSITIVE ENERGY BOUNDARY CONDITION
      IF(ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2) THEN
c
c       get functions at Rmax - dr and Rmax
c
        M=MAX0(LL+1,3)
        T(2)=FK*RMAX
        T(1)=T(2)-FK*DR
        CALL COU(T(1),T(2),ETA,M ,DR,F( 1),F(21),F(41),F(61),F(81))
        SCALE=ABS(F(LL+60))
        T(8)=F(LL+60)/SCALE
        T(7)=F(LL+40)/SCALE
        T(6)=F(LL+20)/SCALE
        T(5)=F(LL   )/SCALE
        D(4)=F(LL+20)
        D(3)=F(LL+60)
c
c     Calculate surface term for width
c
c     calculate Surf = d(k*u'/u)/dk = (u'/u) +kR*[(u''/u) - (u'/u)**2]
c
c      Ratio = (u'/u)
c
      Ratio = (((FL+1.0)**2 + ETA)/t(2) 
     1        - sqrt((FL+1.0)**2 + ETA**2)*F(LL+61)/F(LL+60))/(FL+1.0)
      Surf = ( (1. - T(2)*Ratio)*Ratio 
     1      -T(2)*(1.0 - FL*(FL+1.0)/T(2)**2 - 2.0*ETA/T(2)))/(2.0*fk)
c
c	write(ioutput,'(a,1p4e12.4)')' Surf =',Surf
c
      ELSEIF(ISW.EQ.3.OR.ISW.EQ.4) THEN
        T(8)=0.0
        T(7)=1.0E-15
	SCALE=1.0
      ELSE
        IERR=1
        WRITE(ioutput,'(''0Illegal FISW parameter, ='',i3)')isw
        RETURN
C
      ENDIF
      ELSE
C     NEGATIVE ENERGY BOUNDARY CONDITION
      IF(ISW.EQ.0.OR.ISW.EQ.1) THEN
        T(8)=1.0E-15
        T(7)=T(8)*EXP(FK*DR+ETA* LOG(1.0+DR/RMAX))
	SCALE=1.0
C     BOX BOUNDARY CONDITIONS
      ELSEIF(ISW.EQ.3.OR.ISW.EQ.4) THEN
        T(8)=0.0
        T(7)=1.0E-15
	SCALE=1.0
      ENDIF
      ENDIF
   40 CONTINUE
        width = 0.0
        fnorm0 = 0.5*t(8)**2 + t(7)**2
        F(Kmax  )=T(8)
        F(Kmax-1)=T(7)
        RMAX=RMAX-DR
        K1 = Kmax-2
      turnpt=.true.
      KM=Kmax
C
C     INTEGRATE FROM INFINITY TO TURNING POINT
C
   50 continue
      R=RMAX
      FNCT=0.0
      IF(ISW.NE.2) THEN
      G5=U(2*Kmax-4)-FLP/R**2
      IF(G5.LT.0.0) THEN
        IFLAG=1
      ELSE
        IFLAG=0
      ENDIF
      Q6=1.0+DR2*(U(2*Kmax-2)-FLP/(R+DR)**2)
      Q5=1.0+DR2*G5
      F6=T(8)
      F5=T(7)
      W2=0.0
      FNORM2=fnorm0
C
      DO 100 M=K1,Kmin,-1
      MM = M
      MK=2*M-1
      R=R-DR
      G4=U(MK+1)+V*U(MK  )-FLP/R**2
      Q4=1.0+DR2*G4
      F4=((12.0-10.0*Q5)*F5-Q6*F6)/Q4
      Q6=Q5
      Q5=Q4
      G6=G5
      G5=G4
      if(f6*f5.lt.0.0) fnct=fnct+1.0
      F6=F5
      F5=F4
      F(M)=F4
      IF(G6*G5.LT.0.0) IFLAG=IFLAG+1
      if((.not.turnpt.and.f5.lt.f6).or.IFLAG.GE.2) GO TO 110
      FNORM2=FNORM2+F5**2
      W2=W2+U(MK+2)*F5**2
  100 CONTINUE
C
C     If no classical turning point is found, then find first maximum,
C     then if none, use external matching radius.
C
      if(turnpt) then
        turnpt=.false.
        go to 50
      else
        MM = INT(RM/DR)
      endif
  105 continue
      F5=F(MM  )
      F6=F(MM+1)
  110 CONTINUE
C
C     INTEGRATE FROM ORIGIN TO TURNING POINT
C
      if(isw.lt.2) fnct = 0.0
      KM=MM+1
  120 CONTINUE
      ENDIF
c     End of (Isw ne. 2)
      KS=FL/3.3+2.0
      W1=0.0
      FNORM1=0.0
      F2=0.0
      Q2=0.0
      R =0.0
C
      DO 200 M=Kmin,KM
      MK=M+M-1
      R=R+DR
      Q3=1.0+DR2*(U(MK+1)+V*U(MK  )-FLP/R**2)
      IF(M.GT.KS) THEN
        F3=((12.0-10.0*Q2)*F2-Q1*F1)/Q3
      ELSE
        F3=R**LL
      ENDIF
      Q1=Q2
      Q2=Q3
      F1=F2
      F2=F3
      IF(ICNT.EQ.0) F(M)=F3
      FNORM1=FNORM1+F1**2
      W1=W1+U(MK)*F1**2
C     NODE COUNT
      IF(F1*F2.LT.0.0) FNCT=FNCT+1.0
  200 CONTINUE
      F12=(F1+F2)*0.5
      IF(F1*F2.lt.0.0.and.F12*F2.lt.0.0) FNCT=FNCT-1.0
C
      FN=FNODE-FNCT
      IF(ISW.EQ.2) GO TO 500
      F56=(F5+F6)*0.5
      DET=(F1*F6-F5*F2)/(F12*F56*DR)
      FNORM1=FNORM1/F2**2
      FNORM2=FNORM2/F6**2
      FNORM=FNORM1+FNORM2
      ITER=ITER+1
      IF(ICNT.EQ.0) GO TO 500
      ICNT=ICNT-1
      IF(ICNT.EQ.0) THEN
      WRITE(ioutput,'(''0BOUND STATE SEARCH FAILS IN'',I3,
     1          '' ITERATIONS'')')ICMAX
        IERR=1
      ENDIF
      IF(ICNT.EQ.18) CFACT=CFACT*XFACT
      IF(ICNT.EQ.12) CFACT=CFACT*XFACT
      RM=(FLOAT(KM)-0.5)*DR
      IF((ISW.EQ.3.OR.ISW.EQ.4).AND.E.GT.0.0) THEN
        RSCALE=RMAX
      ELSE
        RSCALE=RM*1.5
      ENDIF
      ESCALE= ABS(((fnode+1.0)*3.0+fl)/(FK2*RSCALE**2))
      DLAMX =DLAMY*ESCALE
        vold=v
        eold=e
      IF(ISW.EQ.0.OR.ISW.EQ.3) THEN
C
C     CHOOSE NEXT GUESS FOR WELL DEPTH
C
        IF(FN.NE.0.0) THEN
c       Node count incorrect
c          dlam= 1.625*cfact*fn*escale
	   dlam= (cfact*((fnode+0.5*fl+0.5)**2
     1                     - (fnct+0.5*fl+0.5)**2) )
          if(dlam.lt.-abs(dv/v)) then
            dlam = abs(dv/v)*sign(one,dlam)
          else
            dlam = min(dlam,dlammax*sign(one,dlam))
          endif
          a='Vn'
        ELSE
          if(det.gt.0.0) then
            dplus =.true.
            if(det.lt.d1) then
              v1=v
              d1=det
            endif
          else
            dminus=.true.
            if(det.gt.d0) then
              v0=v
              d0=det
            endif
          endif
          IF(dplus.and.dminus.and.abs(det).gt.dinterp) THEN
c         Interpolation search
            if(abs(d1/d0).gt.5.0 .or. abs(d1/d0).lt.0.2) then
              DLAM = (v0+v1)/(2.0*v) - 1.0
            else
              DLAM = ((d1*v0-d0*v1)/((d1-d0)*v) - 1.0)
            endif
            a='Vi'
          ELSE
c         Variational search
            DLAM= -DET/(V*DR*(W1/F12**2+W2/F56**2))
            IF(ABS(DLAM).GT.dlammax) DLAM=SIGN(dlammax,DLAM)
            a='Vv'
          ENDIF
        IXNT=IXNT+1
        ENDIF
        flam=1.0+dlam
        V=V*FLAM
        DLAMX=DLAMY
C
C     CHOOSE NEXT GUESS FOR BINDING ENERGY
C
      ELSEIF(ISW.EQ.1.OR.ISW.EQ.4) THEN
        IF(FN.NE.0.0) THEN
c       Node count incorrect
c          DLAM= 1.625*cfact*fn*sign(escale,e)
	   DLAM=  -0.75*cfact*((fnct +0.5*fl+0.5)**2
     1             -(fnode+0.5*fl+0.5)**2)*(3.2/rscale)**2/fk2
          a='En'
        ELSE
          if(det.gt.0.0) then
            dplus =.true.
            if(det.lt.d1) then
              v1=e
              d1=det
            endif
          else
            dminus=.true.
            if(det.gt.d0) then
              v0=e
              d0=det
            endif
          endif
          IF(dplus.and.dminus.and.abs(det).gt.dinterp) THEN
c         Interpolation search
            DLAM= ((d1*v0-d0*v1)/((d1-d0)*e) - 1.0)
            a='Ei'
          ELSE
c         Variational search
            DLAM= -cfact*DET/(DR*FK2*FNORM)
            a='Ev'
          IF(ABS(DLAM).GT.ESCALE) DLAM = SIGN(ESCALE,DLAM)
          ENDIF
        ENDIF
        FLAM=1.0+DLAM
        TEMP=SQRT(ABS(FLAM))
        FK=FK*TEMP
        ETA=ETA/TEMP
        TEMP=FK2*FLAM-FK2
        FK2=FK2+TEMP
        E=E*FLAM
        DO 485 M=1,Kmax
        MK=M+M-1
        U(MK+1)=U(MK+1)+TEMP
  485   CONTINUE
      ENDIF
c     diagnostic printout
c      if(icnt+1.eq.icmax) write(ioutput,8888)
c 8888 format('     Iter',' Nodes',5x,'  Det',9x,'Fnorm',9x,'Old E'
c     1,  9x,'New E',9x,'Old V',9x,'New V',9x,' Flam',9x,'   RM')
c      WRITE(ioutput,8889)A,ITER,FNCT,DET,FNORM,eold,E
c     1  ,vold,V,1.0+DLAM,rm
c 8889 FORMAT(' ',A2,I6,F6.0,8(2X,1PE12.5))
c
      IF(ABS(DET).LT.DETX .AND. ABS(DLAM).LT.DLAMX) THEN
        ICNT=0
      ENDIF
      GO TO 10
C
  500 CONTINUE
C
C     NORMALIZE WAVE FUNCTION
C
      IF(FK2.LT.0.0.OR.ISW.EQ.3.OR.ISW.EQ.4) THEN
C
C     NEGATIVE ENERGY AND BOX B.C.
C
        FNORM=SQRT(FNORM*DR)
      ELSEIF(FK.GT.0.0.AND.(ISW.EQ.0.OR.ISW.EQ.1)) THEN
C     POSITIVE ENERGY AND PI/2 PHASE SHIFT
C
        D(1)=0.0
        D(2)=1.0
        VOL=(FNORM1+FNORM2)*DR*F6**2
c       Volume and surface terms are still divided by SCALE**2
        WIDTH=2.0*E/(FK*(VOL+SURF))
        IF(WIDTH.LT.0.0) then
	  WRITE(ioutput,'(''0 Negative width '',20(1H*)/)')
	endif
        WIDTH=ABS(WIDTH)
        TEMP=( LOG(WIDTH)-2.0* LOG(SCALE))/ELOG10
        I1=TEMP-1.0
        A1=EXP(ELOG10*(TEMP-FLOAT(I1)))
        WRITE(ioutput,9502)A1,I1,VOL*SCALE**2,SURF*SCALE**2
        FNORM=1.0/SQRT(WIDTH*FK/(2.0*E))/F6
C
      ELSEIF(FK.GT.0.0.AND.ISW.EQ.2) THEN
C
C     FIND B.C. FOR E > 0  AND NO SEARCH (ISW = 2)
C
        DET=T(5)*T(8)-T(6)*T(7)
        A1=(F1  *T(8)-F2  *T(7))/DET
        B1=(T(5)*F2  -T(6)*F1  )/DET
        DET=1.0/SQRT(A1**2+B1**2)
        A1=A1*DET
        B1=B1*DET
C     A1=COS(DELTA), B1=SIN(DELTA)
C    (D1,D2)  = (exp(2i*delta) - 1)/(2i) = exp(i*delta)*sin(delta)
        WRITE(ioutput,9501)A1,B1
        D(1)=B1*A1
        D(2)=B1*B1
        FNORM=FK/SCALE
        F2=1.0/DET
        F6=1.0
        DET=0.0
      ENDIF
C
      TEMP=1.0/(F2*FNORM)
      R=0.0
      DO 510 M=1,KM
      R=R+DR
      F(M)=F(M)*TEMP/R
  510 CONTINUE
      IF(KM.LT.Kmax) THEN
        KM=KM+1
        TEMP=1.0/(F6*FNORM)
        DO 520 M=KM,Kmax
        R=R+DR
        F(M)=F(M)*TEMP/R
  520 CONTINUE
      ENDIF
      IF(IPRNT.EQ.0.OR.IERR.NE.0) then
	WRITE(ioutput,9500)V,DET,FNCT,RM,E,ITER
      endif
c
      D(5) =  fk2
      D(6) = (fk2/e)*width/2.0/scale**2
      D(10)=ISW
c
      RETURN
c  
 9500 FORMAT(21X,6HV    =,F9.4,3X,6HDET  =,F9.4,3X,6HNODES=,F9.4,3X,
     1 6HRM   =,F9.4,3X,6HE    =,F9.4,3X,6HITER.=,I4  )
 9501 FORMAT(21X,6HCOSD =,1pe12.4,12x,9H   SIND =,1pe12.4)
 9502 FORMAT(24H0SINGLE PARTICLE WIDTH = ,F7.4,1HE,I3,'  MEV     '
     1,13HVOLUME TERM =  ,1PE12.4,5X,14HSURFACE TERM =  ,1PE12.4)
      END