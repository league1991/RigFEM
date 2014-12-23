% A = loadcrs( 'filename.crs' )
%
% Load a matrix stored under the Compressed Row Storage (CRS) format by
% mat2crs.m 
%

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%
 
function A = loadcrs( filename )

fid = fopen(filename,'r') ;
if fid == -1, error('could not open file'); end;
n  = fscanf(fid,'%u', 1); 
nz = fscanf(fid,'%u', 1);
 
ia = zeros(nz,1); ja = ia; ra = ia;

for i = 1:n+1,
        ia(i) = fscanf(fid,'%u',1);
end;
for i = 1:nz,
        ja(i) = fscanf(fid,'%u',1);
        ra(i) = fscanf(fid,'%f',1);
end;
fclose(fid);

% expand row indices ...
for i = n:-1:1
	for k = 1:ia(i+1)-ia(i)
		ia(nz) = i;
		nz = nz-1;
	end;
end;
A = sparse(ia,ja,ra,n,n); clear ia ja ra;

