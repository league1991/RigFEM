*
*---  sample program illustrating the computation of small matrix
*     exponentials in full with Expokit. Refer to the Expokit
*     documentation for more details about the methods, and 
*     especially the domain of applicability of the Chebyshev scheme.
*
      implicit none
      integer m,i,j,k,ideg,mprint,lda,ldh,lwsp,liwsp,ns,iflag,iexp
      parameter ( ideg=6, lda=50, ldh=lda )
      parameter ( lwsp=4*ldh*ldh+ideg+1, liwsp=ldh )
*
      integer iwsp(liwsp), iseed(4)
      double precision t, A(lda,lda), H(ldh,ldh), y(ldh)
      double precision wsp(lwsp), s1, s2
      complex*16 Hc(ldh,ldh), yc(ldh), wspc(lwsp)
*
      double precision DLARAN
      intrinsic CMPLX, CONJG, MIN, DBLE, IMAG

*-------------------------
*     REAL CASE
*-------------------------

*---  set A = random symmetric negative define matrix ...
      t = 2.0d0
      m = 5
      iseed(1) = 3
      iseed(2) = 7
      iseed(3) = 3
      iseed(4) = 7
      do j = 1,m
         do i = j,m
            A(i,j) = DLARAN( iseed )
            A(j,i) = A(i,j)
         enddo
         A(j,j) = -2.5d0 + A(j,j)
      enddo

*---  maximum number of rows/columns to be printed
      mprint = MIN(5,m)

      print*,"t ="
      print*,t
      print 9000,"REAL SYMMETRIC CASE","*******************************"
      print*,"A = "
      print 9001,( (A(i,j), j=1,mprint), i=1,mprint )
*
 9000 format( /,A,/,A )
 9001 format( <mprint>(1X,D11.4) )
*---  Some compliers (e.g., g77) generate 'Unsupported FORMAT specifier'
*     with the specification above. In this case, simply use this form:
* 9001 format( 5(1X,D11.4) )


*---  Pade ...
      call DGPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

      call DSPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DSPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

*---  Chebyshev
      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DGCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DGCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo

      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DSCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DSCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo

*---  set H = upper Hessenberg part of A ...
      do j = 1,m
         do i = 1,MIN(j+1,m)
            H(i,j) = A(i,j)
         enddo
         do k = i,m
            H(k,j) = 0.0d0
         enddo
      enddo

      print 9000,"REAL UPPER HESSENBERG CASE","************************"
      print*,"H ="
      print 9001,( (H(i,j), j=1,mprint), i=1,mprint )

*---  Pade ...
      call DGPADM( ideg, m, t, H,ldh, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*H) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

*---  Chebyshev
      do i = 1,m
         y(i) = 0.0d0
      enddo
      y(1) = 1.0d0
      call DNCHBV( m,t, A,lda, y, wsp, iflag )
      print 9000,"With DNCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      enddo

*--------------------------------
*     COMPLEX CASE
*--------------------------------

*--   generate the diagonal ...
      do i = 1,m
         s1 = DLARAN( iseed ) - 2.0d0
         Hc(i,i) = s1
      enddo
*---  generate the lower part ...
      do j = 1,m
         do i = j+1,m
            s1 = DLARAN( iseed ) - 0.5d0
            s2 = DLARAN( iseed ) - 0.5d0
            Hc(i,j) = CMPLX( s1,s2 )
         enddo
      enddo
*---  include the conjugate upper part explicitly ...
      do j = 1,m
         do i = j+1,m
            Hc(j,i) = CONJG( Hc(i,j) )
         enddo
      enddo

      print 9000,"COMPLEX HERMITIAN CASE","****************************"
      print 9000," ","Re(H) ="
      print 9001,( (DBLE(Hc(i,j)), j=1,mprint), i=1,mprint )
      print 9000," ","Im(H) ="
      print 9001,( (IMAG(Hc(i,j)), j=1,mprint), i=1,mprint )

      call ZGPADM( ideg, m,t, Hc,ldh, wspc,lwsp,iwsp, iexp, ns,iflag )
      print 9000,"With ZGPADM:","exp(t*H)e_1 ="
      do i = 1,mprint
         print*,wspc(iexp+i-1)
      enddo

      call ZHPADM( ideg, m,t, Hc,ldh, wspc,lwsp,iwsp, iexp, ns,iflag )
      print 9000,"With ZHPADM:","exp(t*H)_e_1 ="
      do i = 1,mprint
         print*,wspc(iexp+i-1)
      enddo

*---  Chebyshev
      do i = 1,m
         yc(i) = 0.0d0
      enddo
      yc(1) = 1.0d0
      call ZGCHBV( m,t, Hc,ldh, yc, wsp, iwsp, iflag )
      print 9000,"With ZGCHBV:","exp(t*H)e_1 ="
      do i = 1,mprint
         print*,yc(i)
      enddo

*---  Note: the Hermitian feature is not useful vis-a-vis Chebyshev

      END




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
