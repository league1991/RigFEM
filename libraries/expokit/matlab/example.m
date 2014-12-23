% MATLAB script illustrating the use of MEXPV.

disp('Loading the matrix ...');
A = loadcrs('../data/c1024.crs');

disp('Computing w = exp(A)e_1 ...');
[n,n] = size(A);
v = eye(n,1);
[w,err] = mexpv(1,A,v);

disp('w(1:10) =');
disp(w(1:10));

disp('err =');
disp(err)
