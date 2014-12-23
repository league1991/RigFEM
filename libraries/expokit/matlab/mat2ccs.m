% mat2ccs( A, 'filename.ccs' )
% 
% Stores the matrix A into a text file under the Compressed Column Storage
% format. The matrix is assumed to be square.
%
% The CCS format holds the matrix column-wise using three column vectors.
% These are usually referred to as JA(1:n+1), IA(1:nz) and A(1:nz).
%
% Upon exit of this script, the first line of the resulting file is 
% composed of n, the order of the matrix, followed by nz, the number of 
% non-zero entries. The rest of the file can be detailed as follows:
% o The first n+1 integers are pointers indicating the beginning 
%   of each column: JA(1:n+1).
% o The remaining part depicts two column vectors:
%   - 	The first column keeps row indices: IA(1:nz).
%   - 	The second column keeps non-zero values: A(1:nz).
%
% Thus, to read the data stored via mat2ccs.m:
% 1. read n, nz
% 2. read JA(1:n+1)
% 3. read (IA(i), A(i), i = 1:nz)

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%

function mat2ccs( A, filename )

%--- retrieve matrix elements ...

[m,n] = size(A);
if ( m~=n ) error('matrix is not square'); end;

[ia,ja,ra] = find(A);
nz = length(ia);
if ( nz==0 ) error('empty matrix'); end;

%--- transform from COO to CCS ...

[ja,jx] = sort(ja);
ia = ia(jx);
ra = ra(jx);

% adjust column pointers
col = zeros(1,n);
for j = 1:nz
  col(ja(j)) = col(ja(j)) + 1;
end;
ptr = 1;
for j = 1:n,
  ja(j) = ptr;
  ptr = ptr + col(j);
end;
ja(n+1) = nz + 1;
ja = ja(1:n+1);

% sort row indices in increasing order
for j = 1:n
  col = ja(j):ja(j+1)-1;
  val = ra(col);
  [ix,jx] = sort( ia(col) );
  ia(col) = ix;
  ra(col) = val(jx);
end

%--- store into the file ...
 
fid = fopen(filename,'w') ;
if ( fid==-1 ) error('could not open file'); end;

% first line
fprintf(fid,'%u  %u\n',[n; nz]);

% column pointers 
fprintf(fid,'%u\n',ja');
 
% remaining part (column indices and values )
fprintf(fid,'%u  %21.15E\n',[ia'; ra']);

fclose(fid);
