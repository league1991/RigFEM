*
*---  sample program illustrating the use of DGPHIV ...
*     Nonhomogeneous problem (Example 6.5 in the Expokit report) ...
*
      implicit none
      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock

*---  matrix data ...
*---  BEWARE: these values must match those in dgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+3)+5*(mmax+3)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm, tmp
      double precision u(nmax),v(nmax),w(nmax), wsave(nmax), wsp(lwsp)

      integer i, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
      double precision DNRM2

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

*---  back to CCS format ...
      call dgcnvr( 'coo','ccs','n', n,n, nz, ia, ja, a, iwsp )

*---  set other input arguments ...
      t = 10.0d0
      tol = 0.0d0
      m = 30
      itrace = 0

*---- First Run: ------------------------------------------------------
***********************************************************************

*---  the operand vector u is set to zero ...
      do i = 1,n
         u(i) = ZERO
      enddo

*---  the operand vector v is set to (1, ..., 1)' ...
      do i = 1,n
         v(i) = ONE
      enddo

*---  compute w = exp(t*A)*v with DGPHIV ...
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         write(*,*) w(i)
      enddo
      call DCOPY( n, w,1, wsave,1 )

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

*----------------------------------------------------------------------|
*----------------------------------------------------------------------|


*---- Second Run ------------------------------------------------------
***********************************************************************

*---  the operand vector u is set to (1, ..., 1)' ...
      do i = 1,n
         u(i) = ONE
      enddo

*---  the operand vector v is set to zero ...
      do i = 1,n
         v(i) = ZERO
      enddo

*---  compute w = t*phi(t*A)*u with DGPHIV ...
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         write(*,*) w(i)
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

*----------------------------------------------------------------------|
*----------------------------------------------------------------------|

*---  now check ||(u + A * w) - wsave||/||wsave||. Due to the specific
*     previous settings, the answer should agree with tol.
*
      call dgccsv( w, v )
      call DAXPY( n, 1.0d0, u,1, v,1 )
      call DAXPY( n, -1.0d0, wsave,1, v,1 )
      tmp = DNRM2( n, v,1 ) / DNRM2( n, wsave, 1 )
      print*
      print*,"relative difference (phi vs. exp) =", tmp
      print*

*---- Third Run -------------------------------------------------------
***********************************************************************

*---  the operand vector v is set to (1, ..., 1)' ...
      do i = 1,n
         v(i) = ONE
      enddo

*---  the operand vector u is set to (1, ..., 1)' ...
      do i = 1,n
         v(i) = ONE
      enddo

*---  compute w = exp(t*A)*v + t*phi(t*A)*u with DGPHIV ...
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
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
*
 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
*----------------------------------------------------------------------|


