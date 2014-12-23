/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _ARPACKSOLVER_H_
#define _ARPACKSOLVER_H_

#include "sparseMatrix.h"

class ARPACKSolver
{
public:
 
  // K * x = lambda * M * x
  // returns the number of converged eigenvalues
  // assumes that both K and M are symmetric, and that M > 0
  // both matrices are given using the entire matrix (not just lower/upper triangle)
  // mode is either "LM" or "SM" (with SM, must also have K > 0)
  // uses mode 2 of ARPACK (regular generalized eigenvalue problem)
  int SolveGenEigReg(SparseMatrix * K, SparseMatrix * M, int numEigenvalues, double * eigenvalues, double * eigenvectors, char * mode = "LM", int numLinearSolverThreads=0, int verbose=1);

  // K * x = lambda * M * x
  // solves for the smallest (in absolute sense) eigenvalues
  // returns the number of converged eigenvalues
  // assumes that both K and M are symmetric, and that M >= 0
  // K can be singular
  // uses mode 3 of ARPACK (shift-inverted generalized eigenvalue problem)
  int SolveGenEigShInv(SparseMatrix * K, SparseMatrix * M, int numEigenvalues, double * eigenvalues, double * eigenvectors, double sigma=0.0, int numLinearSolverThreads=0, int verbose=1);

protected:
};

#endif

