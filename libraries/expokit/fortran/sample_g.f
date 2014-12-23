*
*---  sample program illustrating the use of DGEXPV ...
*     Non-symmetric problem (Example 6.3 in the Expokit report) ...
*
      implicit none
      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock

*---  matrix data ... 
*---  BEWARE: these values should match those in dgmatv.f
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
      double precision t, tol, anorm
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
*
*---  Executable statements ...

*---  load a Harwell-Boeing matrix ...
      n = nmax
      nz = nzmax
      call loadhb( '../data/orani678$', 'coo', n,nz,ia,ja,a, iwsp )

*---  compute the infinite norm of A ...
      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      enddo
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm

*---  the operand vector v is set to (1, ..., 1)^T ...
      do i = 1,n
         v(i) = ONE
      enddo

*---  set other input arguments ...
      t = 10.0d0
      tol = 0.0d0
      m = 30
      itrace = 0

*---  compute w = exp(t*A)v with COO format ...
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcoov, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (COO) has completed:'
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

*--   convert from COO to CCS ...
      call dgcnvr( 'coo','ccs','n', n,n, nz, ia, ja, a, iwsp )

*---  compute w = exp(t*A)v with CCS format ...

      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (CCS) has completed:'
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

*---  convert from CCS to CRS ...
      call dgcnvr( 'ccs','crs','n', n,n, nz, ia, ja, a, iwsp )

*---  compute w = exp(t*A)v with CRS format ...

      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (CRS) has completed:'
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




