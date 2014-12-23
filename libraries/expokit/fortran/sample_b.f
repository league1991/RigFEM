*
*---  sample program illustrating the use of DSEXPV ...
*     Forward-backward problem (Example 6.4 in the Expokit report) ...
*
      implicit none
      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock

*---  matrix data ...
*---- BEWARE: these values must match those in dgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm, tmp
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, j, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
*
*---  Executable statements ...

*---  load the matrix ...
      n = nmax
      nz = nzmax
      call loadhb( '../data/gr3030$', 'crs', n,nz,ia,ja,a, iwsp )

*---  compute the infinite norm of A ...
      anorm = 0.0d0
      do i = 1,n
         tmp = 0.0d0
         do j = ia(i),ia(i+1)-1
            tmp = tmp + ABS( a(j) )
         enddo
         if ( anorm.lt.tmp ) anorm = tmp
      enddo
      write(*,FMT='(A,E8.2)') '||A||_inf= ',anorm

*---  the operand vector v is set to (1, ..., 1)' ...
      do i = 1,n
         v(i) = ONE
      enddo

*---  set other input arguments ...
      t = 1.0d0
      tol = 1.0d-10
      m = 30
      itrace = 0

*---  compute exp(t*A)v with CRS format ...
      tic = clock()
      call DSEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DSEXPV (Forward) has completed:'
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

*---  compute exp(-t*A)w with CRS format ...
      t = -t
      call DCOPY( n, w,1, v,1 )
      tic = clock()
      call DSEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DSEXPV (Backward) has completed:'
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


