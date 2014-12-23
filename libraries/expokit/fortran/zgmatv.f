*-------------------------------NOTE-----------------------------------*
*     This is an accessory to Expokit and it is not intended to be     *
*     complete. It is supplied primarily to ensure an unconstrained    *
*     distribution and portability of the package. The matrix-vector   *
*     multiplication routines supplied here fit the non symmetric      *
*     storage and for a symmetric matrix, the entire (not half) matrix *
*     is required.  If the sparsity pattern is known a priori, it is   *
*     recommended to use the most advantageous format and to devise    *
*     the most advantageous matrix-vector multiplication routine.      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine zgcoov ( x, y )
      implicit none
      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )
 
      do j = 1,n
         y(j) = ZERO
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine zgcrsv ( x, y )
      implicit none
      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Row Storage (CRS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do i = 1,n
         y(i) = ZERO
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine zgccsv( x, y )
      implicit none
      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Column Storage (CCS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do i = 1,n
         y(i) = ZERO
      enddo
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         enddo
      enddo
      end
*----------------------------------------------------------------------|

