c$debug
c***********************************************************************
      PROGRAM DWUCK4
c
c     Main program,  calls ADWCK4 and BDWCK4
c     
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
C
C     STORE STANDARD ANGLE DATA
C
      ANGLE(1)=37.0
      ANGLE(2)=0.0
      ANGLE(3)=5.0
c     Store fundamental constants
      HBARC=197.327
      AMU=931.495
      FSCON=137.036
      CHSQ=HBARC/FSCON
      IBUFF=3200
      IFIRST=0
C
 1001 CONTINUE
C      CALL DW4PC  (IDAT,IFIRST)
C     CALL DW4Vax (IDAT,IFIRST)
      CALL DW4unix(IDAT,IFIRST)
C     CALL DW4ibm (IDAT,IFIRST)
      WRITE(6,9999)IFIRST
c	call second(time)
	write(*,9900)time
	write(6,9900)time
      ifirst = 1
      ibf(3) = 0
C
C     TAPE 2 STORES RADIAL INTEGRALS
      REWIND 2
C
C     TAPE 4 STORES DISTORTED WAVES
      REWIND 4
C

C     READ CARD SET 1

      WRITE(6,9502)
      WRITE(6,9804)IDAT
c
c     icon    Control integers
c     alpha   Alphanumeric run identification
c
      icon(1) = 9
      READ (5,9802,END=2)ICON,ALPHA
      WRITE(6,9803)ICON,ALPHA
    2 continue
      IF(ICON(1).eq.9) then
c
c     close up shop
c
        close (2)
        close (4)
        close (5)
        close (6) 
        STOP ' End of data'
      endif
C
C     DIVISION FOR OVERLAY MAY BE MADE  1-ADWUCK4, 2-BDWUCK4
C
        write(*,'('' '',15a4)')alpha
c        call second(time)
	write(*,'(a)') ' Adwck4 entered'
      CALL ADWUCK4
C     CALL OVERLAY(FNAME,1,0)
      IF(IBF(3).EQ.0) then
c        call second(time)
	write(*,'(a)') ' Bdwck4 entered'
        CALL BDWUCK4
C        CALL OVERLAY(FNAME,2,0)
      endif
      GO TO 1001
c
 9502 FORMAT('0CONTROL INTEGERS')
 9802 FORMAT(20I1,15A4)
 9803 FORMAT(1H ,20I2,4X,15A4)
 9804 FORMAT('  1 2 3 4 5 6 7 8 9 A B C D E F G H I J K'
     1,      '    RUN IDENTIFICATION  ',40X
     2       ,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9900 format(' Elapsed time in seconds =',f7.2)
 9999 FORMAT(I1,'DWUCK4-DISTORTED WAVES U.COLORADO- PC-VERSION'
     1,'  12/Jan /1995' )
      END
