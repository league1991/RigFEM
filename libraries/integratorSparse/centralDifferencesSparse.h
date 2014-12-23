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

A class to timestep dynamics (e.g. nonlinear deformable FEM) 
using the EXPLICIT central differences integrator. Normally, the implicit newmark and
implicit backward Euler classes are recommended for nonlinear deformable FEM
due to stability under large deformations.  See also integratorBase.h.

This implementation follows
WRIGGERS P.: Computational Contact Mechanics. John
Wiley & Sons, Ltd., 2002., page 275

This class supports two damping modes: 
1. Standard Rayleigh model
  D = dampingMassCoef * M + dampingStiffnessCoef * K
     where K is the stiffness matrix AT THE ORIGIN
2. Tangential Rayleigh model (default)
  D = dampingMassCoef * M + dampingStiffnessCoef * K(q)
     where K is the tangential stiffness matrix in the CURRENT deformed configuration q

Mode 1. is computationally faster as it does not require system updates (i.e., 
matrix inversion). Mode 1. is useful, for example, if you want to timestep 
linear modal analysis simulations (stiffness matrix is constant in that case).

Mode 2. gives a better damping model for large deformations, but because the
system matrix changes, requires factoring a linear system anew at each timestep.

In order to use this class, you need to set the timestep very small, or 
else the explicit integrator will go unstable. Roughly speaking, the timestep 
must resolve the highest frequency present in your simulation. 

See also integratorBase.h .

*/

#ifndef _CENTRALDIFFERENCESSPARSE_H_
#define _CENTRALDIFFERENCESSPARSE_H_

#include "integratorBaseSparse.h"
#include "integratorSolverSelection.h"

#ifdef PARDISO
  #include "sparseSolvers.h"
#endif
#ifdef SPOOLES
  #include "sparseSolvers.h"
#endif
#ifdef PCG
  #include "CGSolver.h"
#endif

class CentralDifferencesSparse : public IntegratorBaseSparse
{
public:
  CentralDifferencesSparse(int numDOFs, double timestep, SparseMatrix * massMatrix, ForceModel * forceModel, int numConstrainedDOFs=0, int * constrainedDOFs=NULL, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, int tangentialDampingMode=1, int numSolverThreads=0);

  virtual ~CentralDifferencesSparse();

  inline virtual void SetTimestep(double timestep) { this->timestep = timestep; DecomposeSystemMatrix(); }

  // performs one timestep of simulation
  virtual int DoTimestep(); 

  // sets q, and (optionally) qvel 
  // returns 0 
  virtual int SetState(double * q, double * qvel=NULL);

  // tangentialDampingMode: 
  // 0 = no updates of the damping matrix under deformations (not recommended for large deformations)
  // 1 = update at every timestep (default)
  // k>1 = update every kth timestep
  inline void SetTangentialDampingMode(int tangentialDampingMode) { this->tangentialDampingMode = tangentialDampingMode; }

  virtual void SetInternalForceScalingFactor(double internalForceScalingFactor);

  virtual void ResetToRest();

protected:
  double * rhs;
  double * rhsConstrained;
  SparseMatrix * rayleighDampingMatrix;
  SparseMatrix * tangentStiffnessMatrix;
  SparseMatrix * systemMatrix;
  int tangentialDampingMode;
  int numSolverThreads;
  int timestepIndex;

  void DecomposeSystemMatrix();

  #ifdef PARDISO
    PardisoSolver * pardisoSolver;
  #endif

  #ifdef SPOOLES
    LinearSolver * spoolesSolver;
  #endif

  #ifdef PCG
    CGSolver * jacobiPreconditionedCGSolver;
  #endif
};

#endif

