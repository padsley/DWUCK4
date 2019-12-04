C      ****** Interface routines for a UNIX based system

      SUBROUTINE DW4UNIX(IDAT,ifirst)
      DIMENSION IDAT(6)
C
      if(ifirst.eq.0) then
        OPEN(UNIT=2,STATUS='SCRATCH',FORM='UNFORMATTED')
        OPEN(UNIT=4,STATUS='SCRATCH',FORM='UNFORMATTED')
C
        CALL unixfile(5,6,'DWUCK4 AT YOUR SERVICE')
      endif
C
      CALL unixDAT(IDAT)
C
      RETURN
      END
C                                                  
C     SUBROUTINE UNIXFILE (INPUT, IOUTPUT, TITLE)  
C                                                  
C     @(#)sunfile.f	1.2 90/05/30 10:38:37                 J.J.A. Zalmstra 
C                                                  
C     This subroutine will determine the files associated with STDIN(= unit 5) 
C     and STDOUT (= unit 6) and print a two line title block to STDOUT
C                                                  
C                                                  
      SUBROUTINE UNIXFILE (INPUT, IOUTPUT, TITLE)
      CHARACTER*(*) TITLE                        
C                                                
      integer      getcwd
      character*40 instdat, fdate
      character*40 infile, outfile
      character*40 user
      character*40  cwd
      character*256 arg
      parameter (instdat = 'Mon May 21 11:05:50 1990 ')

      infile = 'Standard Input'
      outfile= 'Standard Output'
      user   = 'unknown'
      cwd    = 'unknown'
c
c Determine the user
c
      call getenv('USER',arg)
      if(lnblnk(arg) .ne. 0) user = arg

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
      if(getcwd(arg) .eq. 0) then
	 call ltrunc(arg, cwd, 40)
      else
	 write(0,*)'Cannot get current directory name'
      endif
      write(ioutput,1030)user,cwd,infile,outfile
      return

      
 99   write(0,1010)arg
        stop
 1010 format('Cannot open file ',a)
 1020 format(1h1,20(1h*),A,20(1h*)//
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

      SUBROUTINE unixDAT(IDAT)                   
      DIMENSION IDAT(6)                          
C                                                
c      CALL IDATE(IDAT(3),IDAT(2),IDAT(1))        
c
c      CALL iTIME(idat(4), idat(5), idat(6))     
c
      RETURN                                    
      END                                       

      SUBROUTINE SECOND(TIME)
C     THIS SUBROUTINE INTERFACES THE SUN SECONDS ROUTINE
C     TO THE SECONDS CALL IN THE PROGRAMS
c     implicit real*8 time
      real  tarray(2), time
      time = etime(tarray)
      RETURN
      END
