%  y = chbv( H, x )
%  CHBV computes the direct action of the matrix exponential on  
%  a vector: y = exp(H)*x. It uses the partial fraction expansion of 
%  the uniform rational Chebyshev approximation of type (14,14). 
%  About 14-digit accuracy is expected if the matrix H is symmetric 
%  negative definite. The algorithm may behave poorly otherwise. 
%
%  See also PADM, EXPOKIT.
 
%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function y = chbv( H, x )
i = sqrt(-1);

% Coefficients and poles of the partial fraction expansion ...
alpha0   =  0.183216998528140087E-11;
alpha(1) =  0.557503973136501826E+02 - i*0.204295038779771857E+03;
alpha(2) = -0.938666838877006739E+02 + i*0.912874896775456363E+02;
alpha(3) =  0.469965415550370835E+02 - i*0.116167609985818103E+02;
alpha(4) = -0.961424200626061065E+01 - i*0.264195613880262669E+01;
alpha(5) =  0.752722063978321642E+00 + i*0.670367365566377770E+00;
alpha(6) = -0.188781253158648576E-01 - i*0.343696176445802414E-01;
alpha(7) =  0.143086431411801849E-03 + i*0.287221133228814096E-03;

theta(1) = -0.562314417475317895E+01 + i*0.119406921611247440E+01;
theta(2) = -0.508934679728216110E+01 + i*0.358882439228376881E+01;
theta(3) = -0.399337136365302569E+01 + i*0.600483209099604664E+01;
theta(4) = -0.226978543095856366E+01 + i*0.846173881758693369E+01;
theta(5) =  0.208756929753827868E+00 + i*0.109912615662209418E+02;
theta(6) =  0.370327340957595652E+01 + i*0.136563731924991884E+02;
theta(7) =  0.889777151877331107E+01 + i*0.166309842834712071E+02;

p = 7;
theta = -theta;
alpha = -alpha;

I = eye(size(H));
y = alpha0*x;
if (isreal(x) & isreal(H)),
  for i = 1:p
    y = y + (H-theta(i)*I) \ (alpha(i)*x);
  end;
  y = real(y);
else
  theta = [theta,conj(theta)];
  alpha = 0.5*[alpha,conj(alpha)];
  for i = 1:2*p
    y = y + (H-theta(i)*I) \ (alpha(i)*x);
  end;
end;

