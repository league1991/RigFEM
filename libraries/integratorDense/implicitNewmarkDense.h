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

A class to timestep dense dynamics (e.g., reduced deformable nonlinear FEM) 
using implicit Newmark. This was the class used for reduced deformable 
nonlinear FEM simulations in references [1,2]. See also newmarkBase.h

Usage note: typical dense system sizes were r=10, r=20, r=30. With larger 
values of r, there will be a computational slowdown.

*/

#ifndef _IMPLICITNEWMARKDENSE_H_
#define _IMPLICITNEWMARKDENSE_H_

#include "integratorBaseDense.h"

class ImplicitNewmarkDense : public virtual IntegratorBaseDense
{
public:

  // there are three choices for the dense solver: positive-definite, symmetric, general
  // positive-definite solver is the most common choice (and cca 2x faster than other options; however, if system matrix becomes singular (rare; only extreme deformations), an error will be issued)
  typedef enum { generalMatrixSolver, symmetricMatrixSolver, positiveDefiniteMatrixSolver } solverType;
  
  ImplicitNewmarkDense(int r, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, solverType solver=positiveDefiniteMatrixSolver, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, int maxIterations = 1, double epsilon = 1E-6, double NewmarkBeta=0.25, double NewmarkGamma=0.5);

  virtual ~ImplicitNewmarkDense();

  inline virtual void SetTimestep(double timestep) { this->timestep = timestep; UpdateAlphas(); }

  // performs one timestep of simulation (returns 0 on success, and 1 on failure)
  // failure can occur, for example, if you are using the positive definite solver and the system matrix has negative eigenvalues
  virtual int DoTimestep(); 

  inline virtual void SetNewmarkBeta(double NewmarkBeta) { this->NewmarkBeta = NewmarkBeta; UpdateAlphas(); }
  inline virtual void SetNewmarkGamma(double NewmarkGamma) { this->NewmarkGamma = NewmarkGamma; UpdateAlphas(); }
  inline void SetMaxIterations(int maxIterations) { this->maxIterations = maxIterations; }
  inline void SetEpsilon(double epsilon) { this->epsilon = epsilon; }

protected:

  ImplicitNewmarkDense(int r, double timestep, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, double NewmarkBeta=0.25, double NewmarkGamma=0.5);

  int symmetricSolver_lwork;
  double * symmetricSolver_work;

  // parameters for implicit Newmark
  double NewmarkBeta, NewmarkGamma;
  double alpha1, alpha2, alpha3, alpha4, alpha5, alpha6;
  double epsilon; 
  int maxIterations;

  solverType solver;

  void UpdateAlphas();
};

#endif

