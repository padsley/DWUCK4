
      SUBROUTINE DW4IBM (IDAT,ifirst)
C
      DIMENSION IDAT(6)
C
      CALL unixDAT(IDAT)
c
      if(ifirst.eq.0) then
        OPEN(UNIT=2,STATUS='SCRATCH',FORM='UNFORMATTED')
        OPEN(UNIT=4,STATUS='SCRATCH',FORM='UNFORMATTED')
C
        CALL unixfile(5,6,'DWUCK4 AT YOUR SERVICE')
      endif
C
      RETURN
      END

C
C     SUBROUTINE UNIXFILE (INPUT, IOUTPUT, TITLE)
C
C     @(#)sunfile.f	1.2 90/05/30 10:38:37               J.J.A. Zalmstra 
C
C     This subroutine will determine the files associated with STDIN(= unit 5) 
C     and STDOUT (= unit 6) and print a two line title block to STDOUT
C
C
      SUBROUTINE UNIXFILE (INPUT, IOUTPUT, TITLE)
C
      CHARACTER*(*) TITLE
C
      integer      getcwd
      character*40 instdat
      character*40 infile, outfile
      character*40 user
      character*40  cwd
      character*256 arg
      character*32 fdate
      parameter (instdat = 'Mon May 21 11:05:50 1990 ')

      infile = 'Standard Input'
      outfile= 'Standard Output'
      user   = 'unknown'
      cwd    = 'unknown'
c
c Determine the user
c
      call getenv('LOGNAME',arg)
      if(lnblnk(arg) .ne. 0) user = arg
      call getenv('PWD    ',arg)
      if(lnblnk(arg) .ne. 0) cwd  = arg

c
c We have read all flags and must now check for input and/or
c outputfilename. 'arg' contains the argument to check, unless
c there are no arguments at all.
c
      nargs = iargc()
      if(nargs .gt. 0) then
	 call getarg(1,arg)
	 if(arg(1:1) .ne. '-') then
	    open(input,file=arg,err=99)
	    call ltrunc(arg, infile, 40)
            rewind input
	 endif
      endif
      if(nargs .gt. 1) then
	 call getarg(2, arg)
	 open(ioutput,file=arg,err=99)
         rewind ioutput
	 call ltrunc(arg, outfile, 40)
      endif
c
c Print title page 
c
      write(6,1020) title,instdat,fdate()
      write(ioutput,1030)user,cwd,infile,outfile
      return

      
 99   write(0,1010)arg
        stop
 1010 format('Cannot open file ',a)
 1020 format(1h ,20(1h*),A,20(1h*)//
     +       ' installed',t20,a,/,' today is ',t20,a,
     +       20(/))
 1030 format(20x,55(1h*)/20x,1h*,t75,1h*/
     +       20x,1h*,' User     : ',a40,t75,1h*/
     +       20x,1h*,' Directory: ',a40,t75,1h*/
     +       20x,1h*,' Input    : ',a40,t75,1h*/
     +       20x,1h*,' Output   : ',a40,t75,1h*/
     +       20x,1h*,t75,1h*,/20x,55(1h*)/1h1)
        end


      subroutine ltrunc(src, dest, maxlen)
c
c copy src to dest but truncate from the left if
c the length of src exceeds maxlen
c
      character*(*) src, dest

      istart = 1
      iend = lnblnk(src)
      if(iend .gt. maxlen) then
	 istart = iend - maxlen - 3
	 dest(1:2) = '<-'
	 dest(3:maxlen) = src(istart:iend)
      else
	 dest = src(istart:iend)
      endif
      return
      end

c***********************************************************************
      SUBROUTINE unixdat(IDAT)
C

	character*32  date, arg
	character*36 month
	integer idat(6)
	common /unixdate/date

	data month/'JanFebMarAprMayJunJulAugSepOctNovDec'/

	iget(i) = ichar(date(i:i))-48

	call system('date > date.dat')
	open(unit = 1, file = 'date.dat',status = 'old')
	read (1, '(a)')date
	close (1)
	call system('rm date.dat')

c	Time 
	idat(4) = iget(12)*10 + iget(13)
	idat(5) = iget(15)*10 + iget(16)
	idat(6) = iget(18)*10 + iget(19)
c	Day
	if(date(09:09) .eq. ' ') date(09:09) = '0' 
	idat(3) = iget(09)*10 + iget(10)

c	Month 
	do 20 i=1,12
	if(date(5:7) .eq. month(3*i-2:3*i)) idat(2) = i
   20	continue
c	Year 
	idat(1) = iget(27)*10 + iget(28)

C
      RETURN
      END

c***********************************************************************
      SUBROUTINE SECOND(TIME)
c
c     This subroutine interfaces IBM 6000 computers  Mclock 
c     routine to the seconds call for elapsed time in the program.
c***********************************************************************
c
c     lcsec is seconds*100
c
	real*8 time
	logical icall
        data icall,lcini/.true., 0/
c
	lcsec = mclock()
	if(icall) then
	  lcini = lcsec
	  icall=.false.
	endif
        time=float(lcsec - lcini)/100.
      return
      end

	function lnblnk(line)
c
	character*80 line
c
	len = 80
	do 10 i=80,1,-1
	if(line(i:i) .ne. ' ') go to 15
	len = i-1
   10	continue
   15	continue
	lnblnk = len
	return
	end

	character*32 function fdate()

	character*32 date
	common/unixdate/date

	fdate = date

	return
	end
