% mat2crs( A, 'filename.crs' )
% 
% Stores the matrix A into a text file under the Compressed Row Storage
% format. The matrix is assumed to be square.
%
% The CRS format holds the matrix row-wise using three column vectors:
% These are usually referred to as IA(1:n+1), JA(1:nz) and A(1:nz).
%
% Upon exit of this script, the first line of the resulting file is 
% composed of n, the order of the matrix, followed by nz, the number of
% non-zero entries. The rest of the file can be detailed as follows:
% o The first n+1 integers are pointers indicating the beginning 
%   of each row: IA(1:n+1).
% o The remaining part depicts two column vectors:
%   - 	The first column keeps column indices: JA(1:nz).
%   - 	The second column keeps non-zero values: A(1:nz).
%
% Thus to read the data stored via mat2crs.m:
% 1. read n, nz
% 2. read IA(1:n+1)
% 3. read (JA(i), A(i), i = 1:nz)

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%

function mat2crs( A, filename )

%--- retrieve matrix elements ...

[m,n] = size(A);
if ( m~=n ) error('matrix is not square'); end;

[ia,ja,ra] = find(A);
nz = length(ia);
if ( nz==0 ) error('empty matrix'); end;

%--- transform from COO to CRS ...

[ia,ix] = sort(ia);
ja = ja(ix);
ra = ra(ix);

% adjust row pointers
row = zeros(1,n);
for i = 1:nz
  row(ia(i)) = row(ia(i)) + 1;
end;
ptr = 1;
for i = 1:n,
  ia(i) = ptr;
  ptr = ptr + row(i);
end;
ia(n+1) = nz + 1;
ia = ia(1:n+1);

% sort column indices in increasing order
for i = 1:n
  row = ia(i):ia(i+1)-1;
  val = ra(row);
  [jx,ix] = sort( ja(row) );
  ja(row) = jx;
  ra(row) = val(ix);
end

%--- store into the file ...
 
fid = fopen(filename,'w') ;
if ( fid==-1 ) error('could not open file'); end;

% first line
fprintf(fid,'%u  %u\n',[n; nz]);

% row pointers 
fprintf(fid,'%u\n',ia');
 
% remaining part (column indices and values )
fprintf(fid,'%u  %21.15E\n',[ja'; ra']);

fclose(fid);
