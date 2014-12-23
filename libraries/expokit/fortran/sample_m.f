*
*---  sample program illustrating the use of DMEXPV and DGEXPV ...
*     Binary Markov Model (Example 6.1 in the Expokit report) ...
*
      implicit none
      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock

*---  matrix data ...
*---  BEWARE: these values must match those in dgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax=5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, j, itrace, iflag
      double precision ZERO, ONE, tmp
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
*
*---  Executable statements ...

*---  load the infinitesimal generator (CRS format)
      open( UNIT=7,STATUS='old',IOSTAT=iflag,FILE='../data/c1024.crs')
      if ( iflag.ne.0 ) stop 'Error - matrix could not be loaded.'
      read( UNIT=7,FMT=* ) n, nz
      if ( nz.gt.nzmax ) stop 'Please increase nzmax.'
      if ( n.gt.nmax ) stop 'Please increase nmax.'
      read( UNIT=7,FMT=* ) (ia(i), i=1,n+1)
      read( UNIT=7,FMT=* ) (ja(i), a(i), i=1,nz)
      close( UNIT=7 )

*---  make sure the infinitesimal generator is transposed,
*     this encoded check prevents from falling in the famous
*     (or rather infamous) `transpose trap' !
*
      do j = 1,2*n
         wsp(j) = 0.0d0
      enddo
      do i = 1,n
         do j = ia(i),ia(i+1)-1
            wsp(i) = wsp(i) + a(j)
            wsp(n+ja(j)) = wsp(n+ja(j)) + a(j)
         enddo
      enddo
      wsp(1) = ABS( wsp(1) )
      wsp(n+1) = ABS( wsp(n+1) )
      do i = 2,n
         wsp(1) = wsp(1) + ABS( wsp(i) )
         wsp(n+1) = wsp(n+1) + ABS( wsp(n+i) )
      enddo
      if ( wsp(n+1).gt.wsp(1) ) then
         print*,'Transposing the input matrix... '
         call tnspos( n, nz, ia, ja, a, iwsp )
      endif

*---  compute the infinite norm of A ...
      anorm = 0.0d0
      do i = 1,n
         tmp = 0.0d0
         do j = ia(i),ia(i+1)-1
            tmp = tmp + ABS( a(j) )
         enddo
         if ( anorm.lt.tmp ) anorm = tmp
      enddo
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm

*---  the operand vector v is set to the first unit basis vector ...
      v(1) = ONE
      do i = 2,n
         v(i) = ZERO
      enddo

*---  set other input arguments ...
      t = 10.0d0
      tol = 1.0d-10
      m = 30
      itrace = 0

*---  compute w = exp(t*A)v with DMEXPV ...
      tic = clock()
      call DMEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DMEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      enddo

*---  display some statistics if desired ...
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

*----------------------------------------------------------------------|
*----------------------------------------------------------------------|

*---  compute w = exp(t*A)v with DGEXPV ...
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      enddo

*---  display some statistics if desired ...
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine tnspos( n, nz, ia, ja, a, iwsp )

      implicit none
      integer          n, nz, ia(nz), ja(nz), iwsp(n)
      double precision a(nz)
*
*-----Purpose----------------------------------------------------------|
*
*---  TNSPOS transposes a CRS matrix. The transposed matrix remains 
*     under the CRS format.
*
*----------------------------------------------------------------------|
*
      integer i, j

      call dgcnvr( 'crs','ccs','n', n,n, nz, ia, ja, a, iwsp )
      do i = 1,nz
         j = ja(i)
         ja(i) = ia(i)
         ia(i) = j
      enddo
      END
*----------------------------------------------------------------------|


