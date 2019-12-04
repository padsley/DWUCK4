
      SUBROUTINE DWPLOT(NTH0,NLOG,SIGPLT,DTHETA,THETA1)
C
C     THIS IS A NEW VERSION OF THE PRINTER PLOT ROUTINE, June 1992 
C
      Parameter (len = 100, pi = 3.141592)
      DIMENSION    SIGPLT(*)
      CHARACTER*1  temp
     1, BLANK /' '/, ASTER /'*'/, FENCE /'|'/, APLUS /'+'/
      CHARACTER*(len+1) dstore 
      CHARACTER*(len+3) astore 
c
	common/ch3files/input, ioutp
C
      ioutput = 6
      if(ioutp .ne. 0) ioutput = ioutp
      nth = abs(nth0)
      if(nlog.gt.0) then
        MAXN= 0 
        DO 30 I=1,NTH
        IF(SIGPLT(I).LE.0.0) SIGPLT(I)=1.0E-20
        SIGPLT(I)= LOG10(SIGPLT(I)) 
        MAXN=MAX0(MAXN,INT(SIGPLT(I) + 100.0)-NLOG)
   30   CONTINUE
        MAXN = MAXN - 100 + 1 
        nlogp= nlog 
      else
        MAXN = -1
        nlogp=  2 
      endif
c     Write header line
        do  40 j=1,len+3
   40   astore(j:j) = BLANK
        do  45 j=0,nlogp
        indx = j*(len/nlogp) + 1
   45   write(astore(indx:indx+2),'(i3)')MAXN +j 
        WRITE(ioutput,9000)astore
c
      DO 100 I=1,NTH
      tint = THETA1 + (i-1)*DTHETA
C     Coose between angle or cos(angle) data
      if(nth0.lt.0) then
        if    (tint.gt. 1.0) tint = 1.0
        if    (tint.lt.-1.0) tint =-1.0
        theta = acos(tint)*180./pi
      else
        theta = tint
      endif
      cs = cos(theta*pi/180.)
      IF(I.EQ.1.OR.I.EQ.NTH) THEN
        TEMP=APLUS
      ELSE
        TEMP=BLANK
      ENDIF
      DSTORE(    1:    1) = APLUS
      DSTORE(len+1:len+1) = APLUS
      DO 70 J=2,LEN
      DSTORE(j:j) = TEMP
   70 CONTINUE
      if(nlog.le.0) dstore(len/2+1:len/2+1) = APLUS
C
      IF(I .EQ. 1 .OR. I .EQ. NTH) THEN
        DO 80 J=0,NLOGP
        indx = j*(len/nlogp) +1
        DSTORE(indx:indx) = FENCE
   80   CONTINUE
      ENDIF
C
      INDX= (SIGPLT(I)-float(MAXN))*float((LEN/NLOGP)) + 1.5
      IF(INDX.GE.1) DSTORE(indx:indx) = ASTER
C
      WRITE(ioutput,9002)THETA,DSTORE,THETA,cs
  100 CONTINUE
C     Write trailer line
        WRITE(ioutput,9000)astore 
c
      RETURN
c
 9000 FORMAT('    THETA ',A103 ,'  THETA ','cos(theta)')
 9002 FORMAT(' ',F8.2 ,'  ',A101,   F8.2,  f10.4 )
      END