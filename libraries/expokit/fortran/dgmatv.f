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
      subroutine dgcoov ( x, y )
      implicit none
      double precision x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j
 
      do j = 1,n
         y(j) = 0.0d0
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine dgcrsv ( x, y )
      implicit none
      double precision x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Row Storage (CRS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j

      do i = 1,n
         y(i) = 0.0d0
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine dgccsv( x, y )
      implicit none
      double precision x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Column Storage (CCS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j

      do i = 1,n
         y(i) = 0.0d0
      enddo
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         enddo
      enddo
      end
*----------------------------------------------------------------------|

