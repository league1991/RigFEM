% mat2coo( A, 'filename.coo' )
% 
% Stores the matrix A into a text file under the sparse COOrdinates storage 
% format. The matrix is assumed to be square.
%
% The COOrdinates storage format holds the matrix by using an UNORDERED set of 
% triplets ( i, j, a_ij ). The number of triplets corresponds to the number
% of non-zero values in the matrix. 
%
% Upon exit of this script, the first line of the resulting file is 
% composed of n, the order of the matrix, followed by nz, the number of
% non-zero entries. The rest of the file depicts three column vectors:
% o 	The first column keeps row indices.
% o 	The second column keeps column indices.
% o 	The third column keeps non-zero values.
%
% Thus to load the data stored via mat2coo.m:
% 1. read n, nz
% 2. read (ia(i), ja(i), a(i), i = 1:nz)

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%

function mat2coo( A, filename )

[m,n] = size(A);
if ( m~=n ) error('matrix is not square'); end;

[ia,ja,ra] = find(A);
nz = length(ia);
if ( nz==0 ) error('empty matrix'); end;
 
fid = fopen(filename,'w') ;
if ( fid==-1 ) error('could not open file'); end;

% first line
fprintf(fid,'%u  %u\n',[n; nz]);
 
% remaining part
fprintf(fid,'%u  %u  %21.15E\n',[ia'; ja'; ra']);

fclose(fid);
