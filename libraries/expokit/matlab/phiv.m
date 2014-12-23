%  [w, err] = phiv( t, A, u, v, tol, m )
%  PHIV computes an approximation of w = exp(t*A)*v + t*phi(t*A)*u
%  for a general matrix A using Krylov subspace projection techniques.
%  Here, phi(z) = (exp(z)-1)/z and w is the solution of the 
%  nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
%  It does not compute the matrix functions in isolation but instead,
%  it computes directly the action of these functions on the 
%  operand vectors. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = phiv( t, A, u, v ) 
%  computes w = exp(t*A)*v + t*phi(t*A)*u using a default tol = 1.0e-7 
%  and m = 30.
%
%  [w, err] = phiv( t, A, u, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = phiv( t, A, u, v, tol )
%  overrides default tolerance.
%
%  [w, err] = phiv( t, A, u, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace.
%
%  Example 1:
%  ----------
%    n = 100;
%    A = rand(n);
%    u = zeros(n,1);
%    v = eye(n,1);
%    w = phiv(1,A,u,v);
%
%
%  Example 2:
%  ----------
%    % generate a random sparse matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%        end;
%    end;
%    u = rand(n,1);
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = phiv(1,A,u,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    x = A\u;
%    w_matlab = expm(full(A))*(v+x)-x;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also EXPV, MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function [w, err] = phiv( t, A, u, v, tol, m )

[n,n] = size(A);
if nargin == 4,
  tol = 1.0e-7;
  m = min(n,30);
end;
if nargin == 5,
  m = min(n,30);
end;

anorm = norm(A,'inf'); 
mxrej = 10;  btol  = 1.0e-7; 
gamma = 0.9; delta = 1.2; 
mb    = m; t_out   = abs(t);
istep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;
sgn = sign(t); istep = 0;
k1 = 3; xm = 1/m; 

w = v;
while t_now < t_out
  V = zeros(n,m+1); H = zeros(m+3,m+3);
  V(:,1) = A*w + u;
  beta = norm(V(:,1));
  V(:,1) = (1/beta)*V(:,1);
  if istep == 0, 
     fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
     t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
     s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s; 
  end;
  istep = istep + 1;
  t_step = min( t_out-t_now,t_new );
  for j = 1:m
     p = A*V(:,j);
     for i = 1:j
        H(i,j) = V(:,i)'*p;
        p = p-H(i,j)*V(:,i);
     end;
     s = norm(p); 
     if s < btol,
        k1 = 0;
        mb = j;
        t_step = t_out-t_now;
        break;
     end;
     H(j+1,j) = s;
     V(:,j+1) = (1/s)*p;
  end; 
  H(1,mb+1) = 1; 
  if k1 ~= 0, 
     H(m+1,m+2) = 1; H(m+2,m+3) = 1;
     h = H(m+1,m); H(m+1,m) = 0;
     avnorm = norm(A*V(:,m+1)); 
  end;
  ireject = 0;
  while ireject <= mxrej,
     mx = mb + max(1,k1);
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
	err_loc = btol; 
        break;
     else
        F(m+1,m+1) = h*F(m,m+2);
        F(m+2,m+1) = h*F(m,m+3);
        p1 = abs( beta*F(m+1,m+1) );
        p2 = abs( beta*F(m+2,m+1) * avnorm );
        if p1 > 10*p2,
           err_loc = p2;
           xm = 1/m;
        elseif p1 > p2,
           err_loc = (p1*p2)/(p1-p2);
           xm = 1/m;
        else
           err_loc = p1;
           xm = 1/(m-1);
        end;
     end;
     if err_loc <= delta * t_step*tol,
        break;
     else
        t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        if ireject == mxrej,
           error('The requested tolerance is too high.');
        end;
        ireject = ireject + 1;
     end;
  end;
  mx = mb + max( 0,k1-2 );
  w = V(:,1:mx)*(beta*F(1:mx,mb+1)) + w;
  
  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
  s = 10^(floor(log10(t_new))-1);
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;

