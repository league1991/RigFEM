%  [w, err, hump] = mexpv( t, A, v, tol, m )
%  MEXPV computes an approximation of w = exp(t*A)*v using Krylov
%  subspace projection techniques. This is a customised version for
%  Markov Chains. This means that a check is done within this code to
%  ensure that the resulting vector w is a probability vector, i.e.,
%  w must have all its components in [0,1], with sum equal to 1.
%  This check is done at some expense and the user may try EXPV
%  which is cheaper since it ignores probability constraints.
%
%  IMPORTANT: The check assumes that the transition rate matrix Q 
%             satisfies Qe = 0, where e = (1,...,1)'. Don't use MEXPV
%             if this condition does not hold. Use EXPV instead.
%             MEXPV/EXPV require A = Q', i.e., the TRANSPOSE of Q. 
%             Failure to remember this leads to wrong results.
%             
%
%  MEXPV does not compute the matrix exponential in isolation but
%  instead, it computes directly the action of the exponential operator
%  on the operand vector. This way of doing so allows for addressing
%  large sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = mexpv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = mexpv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = mexpv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err] = mexpv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1.  However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and
%  ||w(t)||/||v|| are of the same order of magnitude (further details in
%  reference below). Markov chains are usually well-conditioned problems.
%
%  Example:
%  --------
%    % generate a transition rate matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%	 sumj = 0;
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%            sumj = sumj + A(i,j);
%        end;
%	 A(j,j) = A(j,j)-sumj;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
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
%  See also EXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function [w, err, hump] = mexpv( t, A, v, tol, m )

[n,n] = size(A);
if nargin == 3,
  tol = 1.0e-7;
  m = min(n,30);
end;
if nargin == 4,
  m = min(n,30);
end;

anorm = norm(A,'inf'); 
mxrej = 10;  btol  = 1.0e-7; 
gamma = 0.9; delta = 1.2; 
mb    = m; t_out   = abs(t);
istep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s; 
sgn = sign(t); istep = 0;

w = v;
hump = normv;
while t_now < t_out
  istep = istep + 1;
  t_step = min( t_out-t_now,t_new );
  V = zeros(n,m+1); 
  H = zeros(m+2,m+2);

  V(:,1) = (1/beta)*w;
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
  if k1 ~= 0,
     H(m+2,m+1) = 1; 
     avnorm = norm(A*V(:,m+1)); 
  end;
  ireject = 0;
  while ireject <= mxrej,
     mx = mb + k1;
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
	err_loc = btol; 
        break;
     else
        phi1 = abs( beta*F(m+1,1) );
        phi2 = abs( beta*F(m+2,1) * avnorm );
        if phi1 > 10*phi2,
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2,
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
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
  mx = mb + max( 0,k1-1 );
  w = V(:,1:mx)*(beta*F(1:mx,1));       
  beta = norm( w );
  hump = max(hump,beta);

  ineg = 0;
  for i = 1:n
     if w(i) < 0,
        w(i) = 0;
        ineg = ineg + 1;
     end;
  end;
  wnorm = norm(w,1);
  if ineg > 0, w = (1/wnorm)*w; end;
  roundoff = abs(1.0d0-wnorm)/n;
  
  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
  s = 10^(floor(log10(t_new))-1); 
  t_new = ceil(t_new/s) * s;

  err_loc = max(err_loc,roundoff);
  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;