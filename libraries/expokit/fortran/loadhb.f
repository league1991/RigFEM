*----------------------------------------------------------------------|
      subroutine LOADHB( filename, spformat, n, nz, ia, ja, a, iwsp )

      implicit none
      character        filename*80, spformat*3
      integer          n, nz, ia(*), ja(*), iwsp(*)
      double precision a(*)

*---  Purpose ---------------------------------------------------------|
*
*---  LOADHB loads a matrix stored under the Harwell-Boeing format
*     and renders it into the sparse format specified by spformat.
*
*---  Arguments -------------------------------------------------------|
*
*     filename (input)  
*           name of the file containing the matrix.
*           must end with a '$', i.e., filename is in the form: '...$'
*
*     spformat (input)
*           sparse format in which the matrix is forced to be on output
*           accepted values are:
*              'COO' : COOrdinates
*              'CRS' : Compressed Row Storage
*              'CCS' : Compressed Column Storage (default H-B format)
*
*     n (input/output) 
*           On input,  the maximum allowable order
*           On output, the actual order of the matrix loaded
*
*     nz (input/output)
*           On input,  the maximum allowable number of non zero entries
*           On output, the actual number of non zero entries
*
*     ia,ja,a (output)
*           sparse matrix data stored in the format given in spformat 
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly
*           
*     iwsp (workspace) of length >= n
*
*----------------------------------------------------------------------|
*
      character        title*72, key*8, type*3, ptrfmt*16,
     .                 indfmt*16, valfmt*20, rhsfmt*20, rhstyp*1 
      integer          totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow,
     .                 nrhs, nrhsix, i, j, k, io, nn, nnz
      intrinsic        INDEX
*---
      i = INDEX( filename,'$' ) - 1
      if ( i.le.0 ) stop 'in LOADHB. Bad filename'
      open( UNIT=7, STATUS='old', IOstat=io, FILE=filename(1:i) )
      if ( io.ne.0 ) stop 'Could not access Harwell-Boeing matrix'
      read( UNIT=7,FMT=10 ) title, key,
     .     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     .     type, nrow, nn, nnz, nrhs,
     .     ptrfmt, indfmt, valfmt, rhsfmt
      print*, title,'type :',type,' size :',nrow,nn
      print*,'order :',nn,' number of nonzero :',nnz

      if ( nn.gt.n ) stop 'in LOADHB. Please increase n'
      if ( nnz.gt.nz ) stop 'in LOADHB. Please increase nz'

*---  leave if there is no values ... 
      if ( valcrd.le.0 ) stop 'Empty Harwell-Boeing matrix'
      if ( rhscrd.gt.0 ) then
         read( UNIT=7,FMT=11 ) rhstyp, nrhs, nrhsix
         print*,'There is a second hand'
      endif

*---  read data... 
      read( UNIT=7,FMT=ptrfmt ) (ja(i), i = 1,nn+1)
      read( UNIT=7,FMT=indfmt ) (ia(i), i = 1,nnz)
      read( UNIT=7,FMT=valfmt ) (a(i),  i = 1,nnz)
      close( UNIT=7 )

*---  for the sake of experiments, store both parts if symmetric matrix
      if ( type.eq.'RSA' ) then
*---     expand ja indices ...
         k = nnz
         do j = nn,1,-1
            do i = 1,ja(j+1)-ja(j)
               ja(k) = j
               k = k-1
            enddo
         enddo 
*---     insert the other half ...
         k = nnz
         do i = 1,k
            if ( ia(i).ne.ja(i) ) then
               nnz = nnz + 1
               if ( nnz.gt.nz ) stop 'in LOADHB. Please increase nz'
               ia(nnz) = ja(i)
               ja(nnz) = ia(i)
               a(nnz) = a(i)
            endif
         enddo
         type = 'COO'
      else
         type = 'CCS'
      endif
      call dgcnvr( type,spformat,'n', nn,nn, nnz, ia, ja, a, iwsp )
      n = nn
      nz = nnz

      print*,'Harwell-Boeing matrix loaded'
10    format(A72, A8/ 5i14 / A3, 11x, 4i14 / 2a16, 2a20)
11    format(A1, 13x, 2i14)
      end
*----------------------------------------------------------------------|
