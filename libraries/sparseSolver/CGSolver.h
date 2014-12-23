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
  A conjugate gradient solver built on top of the sparse matrix class.
  There are two solver versions: without preconditioning, and with 
  Jacobi preconditioning.

  You can either provide a sparse matrix, or a callback function to
  multiply x |--> A * x .

  The sparse matrix must be symmetric and positive-definite.

  The CG solvers were implemented by following Jonathan Shewchuk's 
  An Introduction to the Conjugate Gradient Method Without the Agonizing Pain:
  http://www.cs.cmu.edu/~jrs/jrspapers.html#cg

  See also sparseMatrix.h .
*/

#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

#include "linearSolver.h"
#include "sparseMatrix.h"

class CGSolver : public LinearSolver
{
public:

  // standard constructor
  CGSolver(SparseMatrix * A);

  // This constructor makes it possible to only provide a
  // "black-box" matrix-vector multiplication routine 
  // (no need to explicitly give the matrix):
  // given x, the routine must compute A * x, and store it into Ax.
  // "data" should not be used/touched by the user-written "black-box" routine.
  // In order to be able to use the Jacobi preconditioner, one can optionally specify the diagonal of the matrix.
  // If diagonal is not specified, SolveLinearSystemWithJacobiPreconditioner will use the identity preconditioner, i.e., it will be identical to SolveLinearSystemWithoutPreconditioner.
  typedef void (*blackBoxProductType)(const void * data, const double * x, double * Ax);
  CGSolver(int n, blackBoxProductType callBackFunction, void * data, double * diagonal=NULL);

  ~CGSolver();

  // solves A * x = b (without preconditioner)
  // A must be symmetric positive-definite
  // input: initial guess (in x)
  // output: solution (in x)
  // "eps" is the convergence criterium: solver converges when the L2 residual errors is less than eps times the initial L2 residual error, must have 0 < eps < 1
  // maximum number of conjugate-gradient iterations is set by "maxIterations"
  // return value is the number of iterations performed 
  // if solver did not converge, the return value will have a negative sign
  int SolveLinearSystemWithoutPreconditioner(double * x, const double * b, double eps=1e-6, int maxIterations=1000, int verbose=0);

  // same as above, except it uses Jacobi preconditioning
  // the employed error metric is M^{-1}-weighted L2 residual error (see Shewchuk)
  int SolveLinearSystemWithJacobiPreconditioner(double * x, const double * b, double eps=1e-6, int maxIterations=1000, int verbose=0);

  virtual int SolveLinearSystem(double * x, const double * b); // implements the virtual method from LinearSolver by calling "SolveLinearSystemWithJacobiPreconditioner" with default parameters

  // computes the dot product of two vectors
  double ComputeDotProduct(double * v1, double * v2); // length of vectors v1, v2 equals numRows (dimension of A)

protected:
  int numRows;
  blackBoxProductType multiplicator;
  void * multiplicatorData;
  SparseMatrix * A; 
  double * r, * d, * q; // terminology from Shewchuk's work
  double * invDiagonal;

  double ComputeTriDotProduct(double * x, double * y, double * z); // sum_i x[i] * y[i] * z[i]
  static void DefaultMultiplicator(const void * data, const double * x, double * Ax);
  void InitBuffers();
};

#endif

