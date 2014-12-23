*----------------------------------------------------------------------|
      subroutine DGCNVR( from,to,diag, nrow,ncol, nz, ia, ja, a, iwsp )
 
      implicit none
      character        from*3, to*3, diag*1
      integer          nrow, ncol, nz, ia(*), ja(*), iwsp(*)
      double precision a(*)

*-----Purpose----------------------------------------------------------|
*
*---  DGCNVR converts a sparse storage format into another sparse
*     storage format. The matrix can be rectangular.
*
*-----Arguments--------------------------------------------------------|
*
*     from   : (input, character*3) the storage format holding the
*              matrix on input. Accepted values are:
*              'COO' : COOrdinates
*              'CRS' : Compressed Row Storage
*              'CCS' : Compressed Column Storage.
*
*     to     : (input, character*3) the storage format holding the
*              matrix on output. Same accepted values as above.
*
*     diag   : (input, character*1) specifies whether or not the 
*              entire diagonal should be lodged, i.e., including 
*              null diagonal entries. (This may be required by
*              some applications that make use of the locations
*              of diagonal elements.)
*              diag = 'N', no specific attention is given to diagonal 
*                          elements.
*              diag = 'D', the entire diagonal is lodged, including 
*                          null elements on the diagonal.
*              if from=to & diag='D', null diagonal entries are
*              explicitly inserted in the actual matrix.
*
*     nrow   : (input) number of rows in the matrix.
*
*     ncol   : (input) number of columns in the matrix.
*
*     nz     : (input/output) number of non-zeros elements.
*              If diag='D' and null diagonal entries are inserted, then 
*              nz is updated on exit and contains the effective number
*              of entries stored. In what follows, nz' (read nz prime)
*              denotes the updated value of nz, nz <= nz' <= nz+n.
*              If diag='N' then nz'=nz.
*
*     ia(*) : (input/output) of declared length .ge. nz'.
*             On input, 
*                if from='CRS', ia(1:nrow+1) contains pointers for the 
*                beginning of each row. 
*                if from='COO', or 'CCS', ia(1:nz) contains row indices.
*             On output,
*                if to='CRS', ia(1:nrow+1) contains pointers for the 
*                beginning of each row. 
*                if to='COO' or 'CCS', ia(1:nz') contains row indices
*                in increasing order in each column.
*
*     ja(*) : (input/output) of declared length .ge. nz'.
*             On input, 
*                if from='CRS', ja(1:ncol+1) contains pointers for the 
*                beginning of each column. 
*                if from='COO' or 'CCS', ja(1:nz) contains col. indices.
*             On output,
*                if to='CRS', ja(1:ncol+1) contains pointers for the 
*                beginning of each column. 
*                if to='COO' or 'CCS', ja(1:nz') contains column indices
*                in increasing order in each row.
*
*     a(*)  : (input/output) On input, a(1:nz) fits the input format
*             and on output, a(1:nz') fits the output format.
*
*   iwsp(*) : (workspace) of declared length .ge. max(nrow,ncol).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     Department of Mathematics, University of Queensland. 
*     Brisbane QLD 4072, Australia. 1996. 
*----------------------------------------------------------------------|

      integer i, j, k, nn, iflag
      character cpy_from*3, cpy_to*3, cpy_diag*1, c*1, S1*3, S2*3 
      logical LSAME

      intrinsic CHAR,ICHAR,LLE,LGE,MIN
      LSAME(S1,S2) = LLE(S1,S2).and.LGE(S1,S2)
*
*---  upper case strings ...
*
      do k = 1,3
         c = from(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_from(k:k) = c
         c = to(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_to(k:k) = c
      enddo
      c = diag(1:1)
      if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
      cpy_diag = c

*---  quick return if possible ...
      if ( LSAME(cpy_from,cpy_to) .and. cpy_diag.eq.'N' ) return
      iflag = 1
      if ( LSAME(cpy_from,'COO') ) iflag = 0
      if ( LSAME(cpy_from,'CRS') ) iflag = 0
      if ( LSAME(cpy_from,'CCS') ) iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in DGCNVR)'
      iflag = 1
      if ( LSAME(cpy_to,'COO') )   iflag = 0
      if ( LSAME(cpy_to,'CRS') )   iflag = 0
      if ( LSAME(cpy_to,'CCS') )   iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in DGCNVR)'

*
*---  transit via COOrdinates format if input is not in COO ...
*
      if ( LSAME(cpy_from,'CRS') ) then
*---     expand ia indices ...
         nn = nz
         do i = nrow,1,-1
            do k = 1,ia(i+1)-ia(i)
               ia(nn) = i
               nn = nn-1
            enddo
         enddo 
      endif
      if ( LSAME(cpy_from,'CCS') ) then
*---     expand ja indices ...
         nn = nz
         do j = ncol,1,-1
            do k = 1,ja(j+1)-ja(j)
               ja(nn) = j
               nn = nn-1
            enddo
         enddo 
      endif
*
*--   if requested, insert diagonal elements even if they are zero...
*
      if ( cpy_diag.eq.'D' ) then
         nn = MIN( nrow,ncol )
         do i = 1,nn
            iwsp(i) = 0
         enddo
         do k = 1,nz
            if ( ia(k).eq.ja(k) ) iwsp(ia(k)) = 1
         enddo
         do i = 1,nn
            if ( iwsp(i).eq.0 ) then
               nz = nz + 1
               ia(nz) = i
               ja(nz) = i
               a(nz) = 0.0d0
            endif
         enddo
      endif

*---  COO convertion ...
      if ( LSAME(cpy_to,'COO') ) return

*---  CRS convertion ...
      if ( LSAME(cpy_to,'CRS') ) call dcmpac( nrow,nz,ia,ja,a, iwsp )

*---  CCS convertion ...
      if ( LSAME(cpy_to,'CCS') ) call dcmpac( ncol,nz,ja,ia,a, iwsp )

      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
*
      subroutine dcmpac( n, nx, ix, ixx, xx, iwsp )
*--   DCMPAC compacts the array ix and sorts ixx and xx
*--   (This is a gateway routine for DGCNVR) ...
*----------------------------------------------------------------------|

      implicit none
      integer          n, nx, ix(nx), ixx(nx), iwsp(n)
      double precision xx(nx)
      integer          k
*
*---  sort ix and carry ixx and xx along ...
*
      call idsrt2( nx, ix, ixx, xx )
*
*---  adjust pointers ...
*
      do k = 1,n
         iwsp(k) = 0
      enddo
      do k = 1,nx
         iwsp(ix(k)) = iwsp(ix(k)) + 1
      enddo
      ix(n+1) = nx + 1
      do k = n,1,-1
         ix(k) = ix(k+1)-iwsp(k)
      enddo
* 
*---  sort ixx in increasing order and carry xx along ...
*
      do k = 1,n
         call idsrt1( iwsp(k), ixx(ix(k)), xx(ix(k)) )
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
*
      subroutine idsrt1( nx, ix, xx )

*---  IDSRT1: indirect sort -- sort ix and carry xx along
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
*----------------------------------------------------------------------|

      implicit none
      integer          nx, ix(nx)
      double precision xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,L
      double precision TX, TTX, R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location 
*---  IT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(I)
         XX(I)   = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*                           
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(J)
         XX(J)   = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            XX(IJ)  = XX(I)
            XX(I)   = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is 
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is 
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         TTX   = XX(L)
         XX(L)  = XX(K)
         XX(K)  = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine idsrt2( nx, ix, ixx, xx )

*---  IDSRT2: indirect sort: sort ix and carry ixx and xx along
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
*----------------------------------------------------------------------|

      implicit none
      integer          nx, ix(nx), ixx(nx)
      double precision xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,JT,JJT,L
      double precision TX, TTX, R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location 
*---  IT, JT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      JT = IXX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(I)
         IXX(I) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(I)
         XX(I)  = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*                           
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(J)
         IXX(J) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(J)
         XX(J)  = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            IXX(IJ)= IXX(I)
            IXX(I) = JT
            JT     = IXX(IJ)
            XX(IJ) = XX(I)
            XX(I)  = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is 
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is 
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         JJT   = IXX(L)
         IXX(L)= IXX(K)
         IXX(K)= JJT
         TTX   = XX(L)
         XX(L) = XX(K)
         XX(K) = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      JT = IXX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      IXX(K+1) = IXX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      IXX(K+1) = JT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
*----------------------------------------------------------------------|
