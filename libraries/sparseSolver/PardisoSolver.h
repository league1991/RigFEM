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

#ifndef _PARDISOSOLVER_H_
#define _PARDISOSOLVER_H_

/*

  Solves A * x = rhs, where A is sparse, usually large, and positive definite.
  Uses a multi-threaded approach to perform the solve.

  The solution is obtained using the the Pardiso library .
  
  Jernej Barbic, MIT, 2007-2009

*/

#include "linearSolver.h"
#include "sparseMatrix.h"

#define MKL_INT int

class PardisoSolver : public LinearSolver
{
public:

  // the constructor will compute the permutation to re-order A 
  // only the topology of A matters for this step
  // A is not modified
  PardisoSolver(const SparseMatrix * A, int numThreads, int positiveDefinite=0, int directIterative=0, int verbose=0);
  virtual ~PardisoSolver();

  MKL_INT ComputeCholeskyDecomposition(const SparseMatrix * A); // perform complete Cholesky factorization

  // solve: A * x = rhs, using the previously computed Cholesky factorization
  // rhs is not modified
  virtual int SolveLinearSystem(double * x, const double * rhs);

  MKL_INT SolveLinearSystemMultipleRHS(double * x, const double * rhs, int numRHS);

  // solve: A * x = rhs, using the direct-iterative solver
  MKL_INT SolveLinearSystemDirectIterative(const SparseMatrix * A, double * x, const double * rhs);

  /*
    For positiveDefinite=0, SolveLinearSystem is equivalent to: ForwardSubstitution + DiagonalSubstitution + BackwardSubstitution
    For positiveDefinite=1, SolveLinearSystem is equivalent to: ForwardSubstitution + BackwardSubstitution (DiagonalSubstitution is not used)
  */

  MKL_INT ForwardSubstitution(double * y, const double * rhs); // L y = rhs
  MKL_INT DiagonalSubstitution(double * v, const double * y);  // D v = y
  MKL_INT BackwardSubstitution(double * x, const double * v);  // L^T x = v

protected:
  int n;
  int numThreads;
  int positiveDefinite;
  int directIterative;
  int verbose;
  double * a;
  int * ia, * ja;

  void *pt[64];
  MKL_INT iparm[64];
  int mtype;
  MKL_INT nrhs; 
  MKL_INT maxfct, mnum, phase, error, msglvl;

  static void DisabledSolverError();
};

#endif

