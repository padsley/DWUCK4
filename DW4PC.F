c***********************************************************************
c   The following four routines should be kept for compilation for a PC
c***********************************************************************
      SUBROUTINE DW4PC(IDAT,ifirst)
c
c     Subroutine for opening files
c
c***********************************************************************
	if(ifirst.eq.0) then
	  open( 2, FORM ='UNFORMATTED', status ='new')
	  open( 4, FORM ='UNFORMATTED', status ='new')
C
	  CALL PCFILE(5,6,'DWUCK4 AT YOUR SERVICE   ')
C
	endif

      CALL PCDATE(IDAT)
      RETURN
      END

c***********************************************************************
      SUBROUTINE PCFILE (INPUT, IOUTPUT, TITLE)
c
c     Opens input and output files
c***********************************************************************
      CHARACTER*25 TITLE,INFILE,OUTFILE
      integer*2 istatus
C
C     open the input and output files
C     -------------------------------
C
	numarg = nargs() - 1
	if(numarg.ge.1) then
	  call getarg(1,infile,istatus)
	else
	  write(*,'(a\)') ' Enter INPUT  file name '
	  read (*,'(a)')  infile
	endif
	if(numarg.ge.2) then
	  call getarg(2,outfile,istatus)
	else
	  write(*,'(a\)') ' Enter OUTPUT file name '
	  read (*,'(a)') outfile
	endif
c
      OPEN ( INPUT,   FILE= INFILE,  STATUS='OLD' )
      OPEN ( IOUTPUT, FILE= OUTFILE, STATUS='NEW' )
C
C
C     output two line title block
C     ---------------------------
C
      WRITE (IOUTPUT,101) TITLE , INFILE 
  101 FORMAT (' ', A25 ,' **** ',' INPUT  FILE NAME = ', A )
C
      WRITE (IOUTPUT,102) OUTFILE
  102 FORMAT ( 33X, 'OUTPUT FILE NAME = ', A ,/)
      RETURN
      END

c***********************************************************************
      SUBROUTINE PCDATE(IDAT)
c
c     time and date routine for MSDOS FORTRAN versions 4.01 and 5.0
c***********************************************************************
c
      integer*2 i1,i2,i3,i4
      DIMENSION IDAT(6)
C                         
c	get date
c
	call getdat(i1,i2,i3)
	idat(1)=i1
	idat(2)=i2
	idat(3)=i3
c
c	get time
c
	call gettim(i1,i2,i3,i4)
	idat(4)=i1
	idat(5)=i2
	idat(6)=i3
c
      RETURN
      END

c***********************************************************************
      SUBROUTINE SECOND(TIME)
c
c     This subroutine interfaces MS Fortran versions 4.0 or 5.0 time 
c     routine to the seconds call for elapsed time in the program.
c***********************************************************************
c
c     lcsec is seconds*100
c
	real*8 time
	integer*2 ihr,imin,isec,ihun, iyr,imon,iday
	integer*4 lcini,lcsec
	logical icall
        data icall,lcini/.true., 0/
c
	call getdat(iyr,imon,iday)
	call gettim(ihr,imin,isec,ihun)
	lcsec=(3600*(24*iday + ihr) + 60*imin + isec)*100 + ihun
	if(icall) then
	  lcini = lcsec
	  icall=.false.
	endif
        time=float(lcsec - lcini)/100.
      return
      end

