/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "integrator" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC     *
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
  A class to timestep large sparse dynamics using implicit Newmark.
  E.g., unreduced nonlinear FEM deformable dynamics.

  See also integratorBase.h .

  This class either uses SPOOLES, PARDISO, or our own Jacobi-preconitioned 
  CG to solve the large sparse linear systems.

  You can switch between these solvers at compile time,
  by modifying the file integratorSolverSelection.h (see below)
  (run-time solver switching would be possible too with more coding).
*/

#ifndef _IMPLICITNEWMARKSPARSE_H_
#define _IMPLICITNEWMARKSPARSE_H_

#ifdef __APPLE__
  #include "TargetConditionals.h"
#endif

// This code supports three different solvers for sparse linear systems of equations:
// SPOOLES, PARDISO, Jacobi-preconditioned Conjugate Gradients
// You must define exactly one of the macros SPOOLES, PARDISO, PCG,
// by changing the file integratorSolverSelection.h .
// PCG is available with our code; look for it in the "sparseMatrix" library (CGSolver.h)
// SPOOLES is available at: http://www.netlib.org/linalg/spooles/spooles.2.2.html
// For PARDISO, the class was tested with the PARDISO implementation from the Intel Math Kernel Library

#include "integratorSolverSelection.h"
#include "sparseMatrix.h"
#include "integratorBaseSparse.h"

#ifdef PARDISO
  #include "sparseSolvers.h"
#endif
#ifdef SPOOLES
  #include "sparseSolvers.h"
#endif
#ifdef PCG
  #include "CGSolver.h"
#endif

class ImplicitNewmarkSparse : public IntegratorBaseSparse
{
public:

  // constrainedDOFs is an integer array of degrees of freedom that are to be fixed to zero (e.g., to permanently fix a vertex in a deformable simulation)
  // constrainedDOFs are 0-indexed (separate DOFs for x,y,z), and must be pre-sorted (ascending)
  // numThreads applies only to the PARDISO solver; if numThreads > 0, the sparse linear solves are multi-threaded; default: 0 (use single-threading)
  ImplicitNewmarkSparse(int r, double timestep, SparseMatrix * massMatrix, ForceModel * forceModel, int positiveDefiniteSolver=0, int numConstrainedDOFs=0, int * constrainedDOFs=NULL, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, int maxIterations = 1, double epsilon = 1E-6, double NewmarkBeta=0.25, double NewmarkGamma=0.5, int numSolverThreads=0); 

  virtual ~ImplicitNewmarkSparse();

  // damping matrix provides damping in addition to mass and stiffness damping (it does not replace it)
  virtual void SetDampingMatrix(SparseMatrix * dampingMatrix);
  inline virtual void SetTimestep(double timestep) { this->timestep = timestep; UpdateAlphas(); }

  // sets q, qvel 
  // automatically computes acceleration assuming zero external force
  // returns 0 on succes, 1 if solver fails to converge
  // note: there are also other state setting routines in the base class
  virtual int SetState(double * q, double * qvel=NULL);

  // performs one step of simulation (returns 0 on sucess, and 1 on failure)
  // failure can occur, for example, if you are using the positive definite solver and the system matrix has negative eigenvalues
  virtual int DoTimestep(); 

  inline void SetNewmarkBeta(double NewmarkBeta) { this->NewmarkBeta = NewmarkBeta; UpdateAlphas(); }
  inline void SetNewmarkGamma(double NewmarkGamma) { this->NewmarkGamma = NewmarkGamma; UpdateAlphas(); }

  // dynamic solver is default (i.e. useStaticSolver=false)
  virtual void UseStaticSolver(bool useStaticSolver);

protected:
  SparseMatrix * rayleighDampingMatrix;
  SparseMatrix * tangentStiffnessMatrix;
  SparseMatrix * systemMatrix;

  double * bufferConstrained;

  // parameters for implicit Newmark
  double NewmarkBeta,NewmarkGamma;
  double alpha1, alpha2, alpha3, alpha4, alpha5, alpha6;
  double epsilon; 
  int maxIterations;

  void UpdateAlphas();
  bool useStaticSolver;

  int positiveDefiniteSolver;
  int numSolverThreads;
  #ifdef PARDISO
    PardisoSolver * pardisoSolver;
  #endif

  #ifdef PCG
    CGSolver * jacobiPreconditionedCGSolver;
  #endif
};

#endif

