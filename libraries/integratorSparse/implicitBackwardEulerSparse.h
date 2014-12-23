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
  A class to timestep large sparse dynamics using implicit backward Euler.
*/

#ifndef _IMPLICITBACKWARDEULERSPARSE_H_
#define _IMPLICITBACKWARDEULERSPARSE_H_

#include "implicitNewmarkSparse.h"

class ImplicitBackwardEulerSparse : public ImplicitNewmarkSparse
{
public:

  // constrainedDOFs is an integer array of degrees of freedom that are to be fixed to zero (e.g., to permanently fix a vertex in a deformable simulation)
  // constrainedDOFs are 0-indexed (separate DOFs for x,y,z), and must be pre-sorted (ascending)
  // numThreads applies only to the PARDISO and SPOOLES solvers; if numThreads > 0, the sparse linear solves are multi-threaded; default: 0 (use single-threading)
  ImplicitBackwardEulerSparse(int r, double timestep, SparseMatrix * massMatrix, ForceModel * forceModel, int positiveDefiniteSolver=0, int numConstrainedDOFs=0, int * constrainedDOFs=NULL, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, int maxIterations = 1, double epsilon = 1E-6, int numSolverThreads=0); 

  virtual ~ImplicitBackwardEulerSparse();

  // sets q, and (optionally) qvel 
  // returns 0 
  virtual int SetState(double * q, double * qvel=NULL);
  virtual int DoTimestep(); 

protected:
};

#endif

