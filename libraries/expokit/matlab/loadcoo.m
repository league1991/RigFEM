% A = loadcoo( 'filename.coo' )
%
% Load a matrix stored under the COOrdinates format by mat2coo.m
%

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%
 
function A = loadcoo( filename )

fid = fopen(filename,'r') ;
if fid == -1, error('could not open file'); end;
n  = fscanf(fid,'%u', 1); 
nz = fscanf(fid,'%u', 1);
 
ia = zeros(nz,1); ja = ia; ra = ia;
for i = 1:nz,
        ia(i) = fscanf(fid,'%u',1);
        ja(i) = fscanf(fid,'%u',1);
        ra(i) = fscanf(fid,'%f',1);
end;
fclose(fid);
 
A = sparse(ia,ja,ra,n,n); clear ia ja ra;

