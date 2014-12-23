/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "sparseSolver" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC   *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*

  Solves A * x = rhs, where A is sparse, usually large, and positive definite.
  Uses a multi-threaded approach to perform the solve.

  The solution is obtained using the the Pardiso library .

  Jernej Barbic, MIT, 2007-2009

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PardisoSolver.h"
#include "sparseSolverAvailability.h"

#ifdef PARDISO_SOLVER_IS_AVAILABLE

// Pardiso solver is available

#ifdef __APPLE__
  #include "TargetConditionals.h"
#endif

#include "mkl.h"

#if defined(_WIN32) || defined(_WIN64)
  #define PARDISO pardiso
#else
  #define PARDISO pardiso
#endif

PardisoSolver::PardisoSolver(const SparseMatrix * A, int numThreads_, int positiveDefinite_, int directIterative_, int verbose_) : numThreads(numThreads_), positiveDefinite(positiveDefinite_), directIterative(directIterative_), verbose(verbose_)
{
  mkl_set_num_threads(numThreads);

  n = A->Getn();

  if (verbose >= 1)
    printf("Converting matrix to Pardiso format...\n");

  int numUpperTriangleEntries = A->GetNumUpperTriangleEntries();

  a = (double*) malloc (sizeof(double) * numUpperTriangleEntries);  
  ia = (int*) malloc (sizeof(int) * (A->GetNumRows() + 1));  
  ja = (int*) malloc (sizeof(int) * numUpperTriangleEntries);  

  int upperTriangleOnly = 1;
  int oneIndexed = 1;
  A->GenerateCompressedRowMajorFormat(a, ia, ja, upperTriangleOnly, oneIndexed);

  if (verbose >= 2)
    printf("numUpperTriEntries: %d\n", numUpperTriangleEntries);

  // permute & do symbolic factorization

  mtype = positiveDefinite ? 2 : -2;
  nrhs = 1; /* Number of right hand sides. */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = verbose >= 1 ? verbose-1 : 0; /* Print statistical information in file */
  error = 0; /* Initialize error flag */

  for (int i = 0; i < 64; i++) 
    iparm[i] = 0;
  iparm[0] = 1; // No solver default
  iparm[1] = 2; // 0=minimum degree ordering, 2=Fill-in reordering from METIS
  iparm[2] = numThreads; // Numbers of processors, value of OMP_NUM_THREADS
  iparm[3] = directIterative ? 62 : 0; //62; // No iterative-direct algorithm
  iparm[4] = 0; // No user fill-in reducing permutation
  iparm[5] = 0; // Write solution into x
  iparm[6] = 0; // Not in use
  iparm[7] = 2; // Max numbers of iterative refinement steps
  iparm[8] = 0; // Not in use
  iparm[9] = 8; //13; // Perturb the pivot elements with 1E-13
  iparm[10] = 0; //1; // Use nonsymmetric permutation and scaling MPS
  iparm[11] = 0; // Not in use
  iparm[12] = 0; // matchings for highly indefinite symmetric matrices
  iparm[13] = 0; // Output: Number of perturbed pivots
  iparm[14] = 0; // Not in use
  iparm[15] = 0; // Not in use
  iparm[16] = 0; // Not in use
  iparm[17] = -1; // Output: Number of nonzeros in the factor LU
  iparm[18] = 0; // no Output: Mflops for LU factorization
  iparm[19] = 0; // Output: Numbers of CG Iterations
  iparm[20] = 1; // pivoting method

/*
  iparm[0] = 0;
  iparm[2] = numThreads;
*/

  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (int i=0; i<64; i++) 
    pt[i] = 0;

  if (verbose >= 1)
    printf("Reordering matrix...\n");

  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, NULL, &nrhs,
          iparm, &msglvl, NULL, NULL, &error);

  if (error != 0)
  {
    printf("Error: Pardiso matrix re-ordering returned non-zero exit code %d.\n", error);
    throw error;
  }

  if (verbose >= 2)
  {
    printf("\nReordering completed ...\n");
    printf("Number of nonzeros in factors = %d\n", iparm[17]);
    printf("Number of factorization MFLOPS = %d\n", iparm[18]);
  }
}

PardisoSolver::~PardisoSolver()
{
  phase = -1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);

  if (error != 0)
    printf("Error: Pardiso Cholesky dealloacation returned non-zero exit code %d.\n", error);

  free(a);
  free(ia);
  free(ja);
}

void PardisoSolver::DisabledSolverError() {}

MKL_INT PardisoSolver::ComputeCholeskyDecomposition(const SparseMatrix * A)
{
  if (directIterative)
    return 0;

  // compute the factorization
  if (verbose >= 1)
    printf("Factoring the %d x %d matrix (%d threads)...\n",n,n, numThreads);

  int upperTriangleOnly = 1;
  int oneIndexed = 1;
  A->GenerateCompressedRowMajorFormat(a, NULL, NULL, upperTriangleOnly, oneIndexed);

  // factor 
  phase = 22;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL,  &error);

  if (error != 0)
    printf("Error: Pardiso Cholesky decomposition returned non-zero exit code %d.\n", error);

  if (verbose >= 1)
    printf("Factorization completed.\n");

  return error;
}

int PardisoSolver::SolveLinearSystem(double * x, const double * rhs)
{
  if (directIterative != 0)
  {
    printf("Error: direct-iterative flag was specified in the constructor (must use SolveLinearSystemDirectIterative routine).\n");
    return 101;
  }

  if (verbose >= 2)
    printf("Solving linear system...(%d threads, using previously computed LU)\n", numThreads);

  phase = 33;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, (double*)rhs, x, &error);

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n"); 

  return (int)error;
}


MKL_INT PardisoSolver::ForwardSubstitution(double * y, const double * rhs)
{
  if (directIterative != 0)
  {
    printf("Error: direct-iterative flag was specified in the constructor (must use SolveLinearSystemDirectIterative routine).\n");
    return 101;
  }

  if (verbose >= 2)
    printf("Performing forward substitution...(%d threads, using previously computed LU)\n", numThreads);

  int maxIterRefinementSteps = iparm[7];
  iparm[7] = 0;

  phase = 331;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, (double*)rhs, y, &error);

  iparm[7] = maxIterRefinementSteps;

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n"); 

  return error;
}

MKL_INT PardisoSolver::DiagonalSubstitution(double * v, const double * y)
{
  if (directIterative != 0)
  {
    printf("Error: direct-iterative flag was specified in the constructor (must use SolveLinearSystemDirectIterative routine).\n");
    return 101;
  }

  if (positiveDefinite == 1)
  {
    printf("Error: diagonal substitution should not be used for positive-definite matrices.\n");
    return 102;
  }

  if (verbose >= 2)
    printf("Performing forward substitution...(%d threads, using previously computed LU)\n", numThreads);

  int maxIterRefinementSteps = iparm[7];
  iparm[7] = 0;

  phase = 332;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, (double*)y, v, &error);

  iparm[7] = maxIterRefinementSteps;

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n"); 

  return error;
}

MKL_INT PardisoSolver::BackwardSubstitution(double * x, const double * y)
{
  if (directIterative != 0)
  {
    printf("Error: direct-iterative flag was specified in the constructor (must use SolveLinearSystemDirectIterative routine).\n");
    return 101;
  }

  if (verbose >= 2)
    printf("Performing forward substitution...(%d threads, using previously computed LU)\n", numThreads);

  int maxIterRefinementSteps = iparm[7];
  iparm[7] = 0;

  phase = 333;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, (double*)y, x, &error);

  iparm[7] = maxIterRefinementSteps;

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n"); 

  return error;
}

MKL_INT PardisoSolver::SolveLinearSystemMultipleRHS(double * x, const double * rhs, int numRHS)
{
  if (directIterative != 0)
  {
    printf("Error: direct-iterative flag was specified in the constructor (must use SolveLinearSystemDirectIterative routine).\n");
    return 101;
  }

  if (verbose >= 2)
    printf("Solving linear system...(%d threads, using previously computed LU)\n", numThreads);

  phase = 33;

  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &numRHS, iparm, &msglvl, (double*)rhs, x, &error);

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n");

  return error;
}

MKL_INT PardisoSolver::SolveLinearSystemDirectIterative(const SparseMatrix * A, double * x, const double * rhs)
{
  if (directIterative != 1)
  {
    printf("Error: direct-iterative flag was not specified in the constructor.\n");
    return 102;
  }

  if (verbose >= 2)
    printf("Solving linear system...(%d threads, direct-iterative)\n", numThreads);

  int upperTriangleOnly = 1;
  int oneIndexed = 1;
  A->GenerateCompressedRowMajorFormat(a, NULL, NULL, upperTriangleOnly, oneIndexed);
  phase = 23;

  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, (double*)rhs, x, &error);

  if (error != 0)
    printf("Error: Pardiso solve returned non-zero exit code %d.\n", error);

  if (verbose >= 2)
    printf("Solve completed.\n"); 

  return error;
}

#else

// Pardiso Solver is not available

PardisoSolver::PardisoSolver(const SparseMatrix * A, int numThreads_, int positiveDefinite_, int directIterative_, int verbose_) : numThreads(numThreads_), positiveDefinite(positiveDefinite_), directIterative(directIterative_), verbose(verbose_) 
{
  DisabledSolverError();
}

PardisoSolver::~PardisoSolver()
{
  DisabledSolverError();
}

void PardisoSolver::DisabledSolverError()
{
  printf("Error: Pardiso solver called, but it has not been installed/compiled/enabled. After installation, enable it in \"sparseSolverAvailability.h\".\n");
  throw 1;
}

MKL_INT PardisoSolver::ComputeCholeskyDecomposition(const SparseMatrix * A)
{
  DisabledSolverError();
  return 1;
}

MKL_INT PardisoSolver::SolveLinearSystem(double * x, const double * rhs)
{
  DisabledSolverError();
  return 1;
}

MKL_INT PardisoSolver::SolveLinearSystemMultipleRHS(double * x, const double * rhs, int numRHS)
{
  DisabledSolverError();
  return 1;
}

MKL_INT PardisoSolver::SolveLinearSystemDirectIterative(const SparseMatrix * A, double * x, const double * rhs)
{
  DisabledSolverError();
  return 1;
}

#endif

