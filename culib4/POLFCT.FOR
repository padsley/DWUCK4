c$debug
      subroutine polfct(max1,maxi,jr,js,theta,Pol,sr
     1          ,iopen20,iout20,nth,ntheta,ALPHA,IDAT)
c     -------------------------------------------------------------
c     max1	initial target     multiplicity
c     maxi	final   target     multiplicity
c     jr	initial projectile multiplicity
c     js	final   projectile multiplicity
c     sr(1...jr, 1...js, 1...max1, 1...maxi)
c     -------------------------------------------------------------
c     Double precision statements  -------------------------------
c      implicit real*8 (a-h,o-z)
c      double complex sr(js,jr,maxi,max1), a(3,3), b(3,3), c(3,3)
c     1   , d(3,3,4), e(3,3,3)
c      real*8 Knn
c     -------------------------------------------------------------
c     Single precision statements  --------------------------------
      complex        sr(js,jr,maxi,max1), a(3,3), b(3,3), c(3,3)
     1   , d(3,3,4), e(3,3,3)
      real*4 Knn
c     -------------------------------------------------------------
      parameter (nsig = 4, npol = 10, nay = 4, nty = 3
     1         , rads = 3.141592/180.)
      logical iopen20, iout20
      
      dimension s(2,2,3), ii(nsig), jj(nsig)
     1         , Pol(npol), Sy(3,3,3), Sij(3,3,4), IDAT(6), ALPHA(15)
     2         ,csig(nsig), Cij(nsig), Dij(nsig), dsig(nsig)
c
c    spin 1/2 matrices for Spin correlation coefficients
c     s    stored as    S_z, S_x and S_y
c
      data s /-1., 0., 0., 1.
     1       , 0., 1., 1., 0.
     2       , 0.,-1., 1., 0./
C
c     SY MATRIX FOR SPIN 0, 1/2 and 1
c
      DATA Sy/0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0
     1,       0.0,-1.0,0.0, 1.0,0.0,0.0, 0.0,0.0,0.0
     2,       0.0,-.707106781,0.0, .707106781,0.0,-.707106781
     3,       0.0,.707106781,0.0/
C
C     SYY = ( 3*SY*SY - S*S  )
C     S22 = (   S^*S^        )*SQRT(3)/4
C     S21 = (   S^*SZ + SZ*S^)*SQRT(3)/2
C     S20 = ( 3*SZ*SZ -  S*S )/SQRT(2)
C
      DATA Sij/-0.5,0.0,-1.5, 0.0,1.0,0.0, -1.5,0.0,-0.5
     1,         0.0,0.0,0.0, 0.0,0.0,0.0, 1.73205081,0.0,0.0
     2,         0.0,0.0,0.0, -1.2247449,0.0,0.0, 0.0,1.2247449,0.0
     3,         0.70710678,0.0,0.0,0.0,-1.4124214,0.0
     4,         0.0,0.0,0.70710678/
C
      data ii/1, 2, 3, 1/
      data jj/1, 2, 3, 2/
      data zero/0.0/
c
c	write(20,'(1p8e12.4)')sr
      cs=cos(theta*rads)
      ss=sin(theta*rads)
      do 20 n=1,nsig
        csig(n)=0.0
        dsig(n)=0.0
   20 continue
      do 30 i=1,npol
	Pol(i) = 0.0
   30 continue
      Dnn = 0.0
      Knn = 0.0
c
      if(jr.gt.3 .or. js.gt.3) go to 1000
c
c     Calculate Dnn = < S_y(initial) * S_y(final) >
c     Calculate pol = < S_y(final  ) >
c     Calculate asy = < S_y(initial) >
c
      do 200 mx=1,max1
      do 190 my=1,maxi
c
      do 180 m =1,jr
      do 170 mp=1,js
      a(mp,m ) = 0.0
      b(mp, m) = 0.0
      c(mp, m) = 0.0
      if(jr .eq. 3) then
      do 115 i=1,nay
      d(mp,m ,i) = 0.0
  115 continue
      endif
      do 130 m1=1,jr
      do 120 m2=1,js
c	Dnn coefficient  -------------------------------------------
      c(mp,m )=c(mp,m ) + sr(m2,m1,my,mx) * cmplx(zero, Sy(m2,mp,js))
     1     *cmplx(zero, Sy(m ,m1,jr))
  120 continue
c	Asymmetry        --------------------------------------------
      b(mp,m )=b(mp,m ) + sr(mp,m1,my,mx) * cmplx(zero, Sy(m ,m1,jr))
      if(jr .eq. 3) then
      do 125 i=1,nay
      d(mp,m ,i)=d(mp,m ,i) + sr(mp,m1,my,mx) *         Sij(m ,m1,i) 
  125 continue
      endif
  130 continue
      do 140 m2=1,js
c	Polarization     --------------------------------------------
      a(mp,m )=a(mp,m ) + sr(m2,m ,my,mx) * cmplx(zero, Sy(m2,mp,js))
      if(js .eq. 3) then
      do 135 i=1,nty
      e(mp,m ,i)=e(mp,m ,i) + sr(m2,m ,my,mx) *         Sij(m2,mp,1) 
  135 continue
      endif
  140 continue
c
      Pol(1)  =Pol(1) + conjg(sr(mp,m ,my,mx)) * sr(mp,m ,my,mx)
      Pol(2)  =Pol(2) + conjg(sr(mp,m ,my,mx)) * a(mp,m )
      Pol(3)  =Pol(3) + conjg(sr(mp,m ,my,mx)) * b(mp,m )
      if(jr .eq. 3) then
      do 160 i=1,nay
      pol(i+3)=Pol(i+3) + conjg(sr(mp,m ,my,mx)) * d(mp,m ,i)
  160 continue
      endif
      if(js .eq. 3) then
      do 165 i=1,nty
      pol(i+7)=Pol(i+7) + conjg(sr(mp,m ,my,mx)) * e(mp,m ,i)
  165 continue
      endif
      Dnn   =Dnn    + conjg(sr(mp,m ,my,mx)) * c(mp,m )
  170 continue
  180 continue
  190 continue
  200 continue
c
      if(Pol(1) .eq. 0.0) go to 1000
      IF(iout20) THEN
      if(jr.eq.2 .and. maxi.eq.2) then
c
c     Calculate Knn = < S_y(initial) * I_y(final) >
c
      do 300 mx=1,max1
      do 290 mp=1,js
c
      do 280 my=1,maxi
      do 270 m = 1,jr
      c(m ,my) = 0.0
      do 260 m1=1,maxi
      do 250 m2=1,jr
c     Knn coefficient
      c(m ,my)=c(m ,my) + sr(mp,m1,m2,mx) * cmplx(zero, Sy(my,m2,2))
     1     *cmplx(zero, Sy(m1,m ,2))
  250 continue
  260 continue
      Knn = Knn + conjg(sr(mp,m ,my,mx)) * c(m ,my)
  270 continue
  280 continue
  290 continue
  300 continue
      endif
c
      if(jr .eq. 2 .and. js .eq.2) then
c
c     Calculate Dij = < S_i(initial) * S_j(final) >
c

      do 600 mx=1,max1
      do 580 my=1,maxi
      do 500 n=1,nsig
        i1=ii(n)
        j1=jj(n)
      do 490 m =1,jr
      do 480 mp=1,js
        a(mp,m )=0.0
      do 440 m1=1,jr
      do 420 m2=1,js
      a(mp,m )=a(mp,m ) + sr(m2,m1,my,mx) * s(m2,mp,i1)*s(m ,m1,j1)
  420 continue
  440 continue
c
c     Dij correlation coefficients
      dsig(n)=dsig(n)+conjg(sr(mp,m ,my,mx))*a(mp,m )
  480 continue
  490 continue
  500 continue
  580 continue
  600 continue
      endif
c
      Dsig(3) = -Dsig(3)
      do 610 n=1,nsig
      dsig(n) = dsig(n)/Pol(1)
  610 continue
c
      if(js .eq. 2 .and. maxi .eq. 2) then
c
c     Spin correlation coefficients (final state target and projectile)
c     Calculate Cij = < S_y(final) * I_y(final) >
c
      do 800 mx=1,max1
      do 780 m =1,jr
c
      do 700 n=1,nsig
        i1=ii(n)
        j1=jj(n)
      do 690 my=1,maxi
      do 680 mp=1,js
        a(mp,my)=0.0
      do 640 m1=1,maxi
      do 620 m2=1,js
      a(mp,my)=a(mp,my)+sr(m2,m ,m1,mx)*s(m2,mp,i1)*s(m1,my,j1)
  620 continue
  640 continue
c
      csig(n)=csig(n)+conjg(sr(mp,m ,my,mx))*a(mp,my)
  680 continue
  690 continue
  700 continue
c
  780 continue
  800 continue
      do 820 n=1,nsig
	csig(n) = csig(n)/Pol(1)
  820 continue
        csig(3)=-csig(3)
c
c     rotate operators to outgoing particle direction
c
c	Minus signs on C_zz, C_xx and C_xz make output agree with the data
c	where the z and x axes for the target are in opposite directions
c	to those of the projectile
c
      Cij(3) = -(csig(1)*cs**2 + csig(2)*ss**2 + 2.0*csig(4)*cs*ss)
      Cij(1) = -(csig(1)*ss**2 + csig(2)*cs**2 - 2.0*csig(4)*cs*ss)
      Cij(4) = -(csig(4)*(cs**2-ss**2) + (csig(2)-csig(1))*cs*ss)
      Cij(2) =   csig(3)
      endif
c
c     Singlet fraction
      ssum = (1.0-(csig(1)+csig(2)+csig(3)))/4.0
      Dij(3) = dsig(1)
      Dij(1) = dsig(2)
      Dij(4) = dsig(4)
      Dij(2) = dsig(3)
      ENDIF
  900 continue
c
      do 980 i=2,npol
        Pol(i)  = Pol(i)/Pol(1)
  980 continue
        Dnn     = Dnn   /Pol(1)
        Knn     = Knn   /Pol(1)
        Pol(1)  = Pol(1)/float(max1*jr)
 1000 continue
c
c     --------------------------------------------------------
c       output to disk file 20 and file 21
c
	if(iopen20) then
	  open(unit = 20, file = 'for020.dat', status = 'unknown')
	  open(unit = 21, file = 'for021.dat', status = 'unknown')
	  iopen20 = .false.
	endif
c
        if(iout20) then
c       Write header to file
        if(nth .eq. 1) then
          WRITE(20,9010)ALPHA,IDAT
          write(20,9020) ntheta
          WRITE(21,9010)ALPHA,IDAT
          write(21,9021) ntheta
        endif
c
          write(20,'(2(0pf8.3,1h,), 1pe12.4, 9(1h,,0pf8.4))')
     1    theta, cs, (Pol(i),i=1,3), Dnn, Knn
c
          write(21,'(0pf8.3,1h,,0pf8.3,      9(1h,,0pf8.4))')
     1    theta, cs, Cij, Dij, ssum
        endif
c
	return
c 
 9010 FORMAT(' (',15A4,I4.2,2(1H/,I2.2),I4.2,2(1H:,I2.2))
 9020 FORMAT(' (',i2,',angle   cos[th]    Sigma       Pol      Asy'
     1    ,'      Dnn      Knn')
 9021 FORMAT(' (',i2,',angle cos[th]     Cxx      Cyy      Czz'
     1    ,'      Cxz      Dxx      Dyy      Dzz      Dxz   fsingl')
c
      end