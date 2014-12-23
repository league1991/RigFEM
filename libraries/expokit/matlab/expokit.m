%  EXPOKIT is a set of user-friendly routines (in FORTRAN 77 and MATLAB)
%  aimed at computing matrix exponentials. More precisely, it computes
%  either a small matrix exponential in full, the action of a large
%  sparse matrix exponential on an operand vector, or the solution of a
%  system of linear ODEs with constant inhomogeneity. The backbone of
%  the sparse routines consists of Krylov subspace projection methods
%  (Arnoldi and Lanczos processes) and that is why the toolkit is capable
%  of coping with sparse matrices of very large dimension. The software
%  handles real and complex matrices and provides specific routines for
%  symmetric and Hermitian matrices. When dealing with Markov chains,
%  the computation of the matrix exponential is subject to probabilistic
%  constraints.  In addition to addressing general matrix exponentials,
%  a distinct attention is assigned to the computation of transient
%  states of Markov chains.
%  
%  See http://www.maths.uq.edu.au/expokit
%  
%  Notes:
%  ======
%  1. w(t) = exp(t*A)v is the analytic solution of the homogeneous 
%     ODE problem:
%                 w' = Aw,  w(0) = v
%  
%  2. w(t) = exp(t*A)*v + t*phi(t*A)u, where phi(z) = (exp(z)-1)/z,
%     is the analytic solution of the nonhomogeneous ODE problem:
%                 w' = Aw + u,  w(0) = v
%  
%  
%  The MATLAB version implements the following subset:
%  
%  --- computational scripts:
%  
%  chbv.m      Chebyshev algorithm for w = exp(t*A)*v
%  padm.m      irreducible Pade algorithm for exp(t*A)
%  phiv.m      Krylov algorithm for w = exp(t*A)*v + t*phi(t*A)*u
%  expv.m      Krylov algorithm for w = exp(t*A)*v
%  mexpv.m     Markov algorithm for w = exp(t*A)*v
%  
%  
%  --- utilities:
%  
%  mat2ccs.m  save matrix into an ascii file under CCS format
%  mat2coo.m  save matrix into an ascii file under COO format
%  mat2crs.m  save matrix into an ascii file under CRS format
%  	   
%  loadccs.m  load matrix from a file under CCS format
%  loadcoo.m  load matrix from a file under COO format
%  loadcrs.m  load matrix from a file under CRS format
%  	 
%  example.m  simple illustrative example
%  
%  Reference:
%  ==========
%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998


