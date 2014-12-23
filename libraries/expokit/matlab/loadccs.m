% A = loadccs( 'filename.ccs' )
%
% Load a matrix stored under the Compressed Column Storage (CCS) format by 
% mat2ccs.m
%

% Roger B. Sidje (rbs@maths.uq.edu.au), 1996.
%
 
function A = loadccs( filename )

fid = fopen(filename,'r') ;
if fid == -1, error('could not open file'); end;
n  = fscanf(fid,'%u', 1); 
nz = fscanf(fid,'%u', 1);
 
ia = zeros(nz,1); ja = ia; ra = ia;

for j = 1:n+1,
        ja(j) = fscanf(fid,'%u',1);
end;
for i = 1:nz,
        ia(i) = fscanf(fid,'%u',1);
        ra(i) = fscanf(fid,'%f',1);
end;
fclose(fid);

% expand column indices ...
for j = n:-1:1
	for k = 1:ja(j+1)-ja(j)
		ja(nz) = j;
		nz = nz-1;
	end;
end;
A = sparse(ia,ja,ra,n,n); clear ia ja ra;


