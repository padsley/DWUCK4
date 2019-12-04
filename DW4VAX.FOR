c      ****** Interface routine for a VAX with VMS

c***********************************************************************
c   The following five routines should be kept for compilation for a VAX
c   Choose the version of VAXFILE for either input/output routine names
c   to be entered from a command file or from within the program.
c***********************************************************************
      SUBROUTINE DW4VAX(IDAT,ifirst)
c
c     opens files - input, output and scratch files 2 and 4
c***********************************************************************
c
      DIMENSION IDAT(6)
C
      if(ifirst.eq.0) then
       open(unit=2,status='scratch',FORM='UNFORMATTED')
       open(unit=4,status='scratch',FORM='UNFORMATTED')
C
      CALL VAXFILE(5,6,'DWUCK4 AT YOUR SERVICE')
C
      endif
	call vaxdat(idat)
      RETURN
      END

c***********************************************************************
C
C     SUBROUTINE VAXFILE (INPUT, IOUTPUT, TITLE)
C
C                                        Dale E. Prull  03/02/83
C
C     This subroutine will open the following files from a command file:
C
C                    'INFILE'  -- logical unit INPUT
C                    'OUTFILE' -- logical unit IOUTPUT
C
C     The logical name translation for these files is obtained, and
C     a two line "title block" is written to 'OUTFILE' containing:
C
C                    TITLE , date, and time
C                    'INFILE', and 'OUTFILE' logical name translations
C
C
c***********************************************************************
      SUBROUTINE VAXFILE (INPUT, IOUTPUT, TITLE)
c
c     Opens input and output files from within the program
c***********************************************************************
C
      CHARACTER*(*) TITLE
C
      CHARACTER*9 DATEBUF
      CHARACTER*8 TIMEBUF
      CHARACTER*63 NAMEIN,NAMEOUT
      INTEGER*2 LENIN,LENOUT
      INTEGER*4 SYS$TRNLOG
C
C     open the input and output files
C     -------------------------------
C
      OPEN ( UNIT=INPUT,   FILE='INFILE',  STATUS='OLD', READONLY )
      OPEN ( UNIT=IOUTPUT, FILE='OUTFILE', STATUS='NEW',
     $                     CARRIAGECONTROL='FORTRAN' )
C
C     find the logical name translation
C     ---------------------------------
C
      ISTAT=SYS$TRNLOG ( 'INFILE',  LENIN,  NAMEIN ,,, )
      ISTAT=SYS$TRNLOG ( 'OUTFILE', LENOUT, NAMEOUT,,, )
C
C     get the date and time
C     ---------------------
C
      CALL DATE ( DATEBUF )
      CALL TIME ( TIMEBUF )
C
C     output two line title block
C     ---------------------------
C
      WRITE (IOUTPUT,101) TITLE , NAMEIN(1:LENIN) , DATEBUF, TIMEBUF
  101 FORMAT (' ', A25 ,' **** ',' INPUT  FILE NAME = ', A ,
     + T105, A, 2X, A )
C
      WRITE (IOUTPUT,102) NAMEOUT(1:LENOUT)
  102 FORMAT ( 33X, 'OUTPUT FILE NAME = ', A ,/)
      RETURN
      END

c***********************************************************************
      SUBROUTINE VAXDAT(IDAT)
c***********************************************************************
      character*1 STR(8)
      character*8 STX
      DIMENSION IDAT(6)
      equivalence (str,stx)
C
      CALL IDATE(IDAT(3),IDAT(2),IDAT(1))
C
      CALL TIME(STX)
      DO 160 I=1,3
      J=3*I-2
      IDAT(I+3)=(ICHAR(STR(J))-48)*10+ICHAR(STR(J+1))-48
  160 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE SECOND(TIME)
c
C     THIS SUBROUTINE INTERFACES THE VAX SECONDS ROUTINE
C     TO THE SECONDS CALL IN THE PROGRAMS
c***********************************************************************
      REAL*8 TIME
      DATA ICALL/0/
      IF(ICALL.EQ.0) CALL LIB$INIT_TIMER()
        ICALL=1
        CALL LIB$STAT_TIMER(2,LCSEC)
        TIME=0.01*FLOAT(LCSEC)
      RETURN
      END


