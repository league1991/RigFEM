*
*---  sample program illustrating the use of ZGEXPV and ZHEXPV ...
*     Hermitian problem (Example 6.2 in the Expokit report) ...
*
      implicit none
      external zgcoov, zgcrsv, zgccsv

      double precision tic, tac, clock

*---  matrix data ...
*---  BEWARE: these values must match those in zgmatv.f
      integer n, nz, nmax, nzmax
      parameter( nmax = 5500, nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = mmax+2 )
      integer iwsp(liwsp)
      double precision t, tol, anorm, s1, s2
      complex*16 v(nmax), w(nmax), wsp(lwsp)

      integer i, j, nnz, itrace, iflag, iseed(4)
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

      double precision DLARAN
      intrinsic ABS, CMPLX, CONJG, DBLE
*
*---  Executable statements ...

*---  load a symmetric pattern ...
      n = nmax
      nz = nzmax/2
      call getpat( '../data/bcspwr10$', n, nz, ia, ja )

*---  for the purpose of the experiments, expand to COOrdinates ...
      nnz = nz
      do j = n,1,-1
         do i = 1,ja(j+1)-ja(j)
            ja(nnz) = j
            nnz = nnz-1
         enddo
      enddo
*---  fill-in an Hermitian matrix -- the conjugate part is included
      iseed(1) = 0
      iseed(2) = 0
      iseed(3) = 0
      iseed(4) = 5
      nnz = nz
      do i = 1,nz
         if ( ia(i).ne.ja(i) ) then
            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
            s2 = 10.0d0*DLARAN( iseed ) - 5.0d0
            a(i) = CMPLX( s1,s2 )
            nnz = nnz + 1
            a(nnz) = CONJG( a(i) )
            ia(nnz) = ja(i)
            ja(nnz) = ia(i)
        else
            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
            a(i) = CMPLX( s1,0.0d0 )
         endif
      enddo
      nz = nnz

*---  compute the infinite norm of A ...
      do i = 1,n
         wsp(i) = ZERO
      enddo
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      enddo
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.DBLE(wsp(i)) ) anorm =  wsp(i)
      enddo

**--   convert from COO to CRS ...
*      call zgcnvr( 'coo','crs','n', n,n, nz, ia, ja, a, iwsp )
*
**---  compute the infinite norm of A ...
*      anorm = 0.0d0
*      do i = 1,n
*         s1 = 0.0d0
*         do j = ia(i),ia(i+1)-1
*            s1 = s1 + ABS( a(j) )
*         enddo
*         if ( anorm.lt.tmp ) anorm = s1
*      enddo

      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm

*---  the operand vector v is set to e_1 + e_n ...
      do i = 1,n
         v(i) = ZERO
      enddo
      v(1) = ONE
      v(n) = ONE

*---  set other input arguments ...
      t = 1.0d0
      tol = 1.0d-5
      m = 30
      itrace = 0

*---  compute w = exp(t*A)v with ZGEXPV ...
      tic = clock()
      call ZGEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'ZGEXPV has completed:'
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
      print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
      print 9002,'step_min  = ',DBLE( wsp(1) )
      print 9002,'step_max  = ',DBLE( wsp(2) )
      print 9002,'max_round = ',DBLE( wsp(3) )
      print 9002,'sum_round = ',DBLE( wsp(4) )
      print 9002,'max_error = ',DBLE( wsp(5) )
      print 9002,'sum_error = ',DBLE( wsp(6) )
      print 9002,'hump      = ',DBLE( wsp(9) )
      print 9002,'scale-norm= ',DBLE( wsp(10) )

*----------------------------------------------------------------------|
*----------------------------------------------------------------------|

*---  compute w = exp(t*A)v with ZHEXPV ...
      tic = clock()
      call ZHEXPV( n, m, t,v,w, tol, anorm,
     .             wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      tac = clock()

*---
      print 9001,'----------------------------------------------------'
      print 9001,'ZHEXPV has completed:'
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
      print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
      print 9002,'step_min  = ',DBLE( wsp(1) )
      print 9002,'step_max  = ',DBLE( wsp(2) )
      print 9002,'max_round = ',DBLE( wsp(3) )
      print 9002,'sum_round = ',DBLE( wsp(4) )
      print 9002,'max_error = ',DBLE( wsp(5) )
      print 9002,'sum_error = ',DBLE( wsp(6) )
      print 9002,'hump      = ',DBLE( wsp(9) )
      print 9002,'scale-norm= ',DBLE( wsp(10) )

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine getpat( filename, n, nz, ia, ja )
*---  load a Harwell-Boeing pattern ...

*---  arguments are fully described in LOADHB
      implicit none
      character        filename*80
      integer          n, nz
      integer          ia(*), ja(*)
*---
      character        title*72, key*8, type*3, ptrfmt*16,
     .                 indfmt*16, valfmt*20, rhsfmt*20, rhstyp*1
      integer          totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow,
     .                 nrhs, nrhsix, i, io, nn, nnz
*---
      i = INDEX( filename,'$' ) - 1
      if ( i.le.0 ) stop 'in GETPAT. Bad filename'
      open( UNIT=7, STATUS='old', IOstat=io, FILE=filename(1:i) )
      if ( io.ne.0 ) stop 'Could not access Harwell-Boeing matrix'
      read( UNIT=7,FMT=10 ) title, key,
     .     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     .     type, nrow, nn, nnz, nrhs,
     .     ptrfmt, indfmt, valfmt, rhsfmt
      print*, title,'type :',type,' size :',nrow,nn
      print*,'order :',nn,' number of nonzero :',nnz

      if ( nn.gt.n ) stop 'Please increase nmax.'
      if ( nnz.gt.nz ) stop 'Please increase nzmax.'
      n = nn
      nz = nnz
      if ( rhscrd.gt.0 ) then
         read( UNIT=7,FMT=11 ) rhstyp, nrhs, nrhsix
         print*,'There is a second hand'
      endif

*---  read data... 
      read( UNIT=7,FMT=ptrfmt ) (ja(i), i = 1,n+1)
      read( UNIT=7,FMT=indfmt ) (ia(i), i = 1,nz)

      close( UNIT=7 )
      print*,'Harwell-Boeing pattern loaded'
10    format(A72, A8/ 5i14 / A3, 11x, 4i14 / 2a16, 2a20)
11    format(A1, 13x, 2i14)
      end
*----------------------------------------------------------------------|
      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      DLARAN = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $         ( DBLE( IT4 ) ) ) ) )
      RETURN
*
*     End of DLARAN
*
      END
*----------------------------------------------------------------------|
