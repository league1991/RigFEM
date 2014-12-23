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

A class to timestep dense dynamics (e.g. reduced nonlinear deformable FEM) 
using explicit central differences. Normally, the ImplicitNewmark class is 
recommended for reduced nonlinear deformable FEM. See also newmarkBase.h.

This class supports two damping modes: 
1. Standard Rayleigh model
  D = dampingMassCoef * M + dampingStiffnessCoef * K
     where K is the stiffness matrix AT THE ORIGIN
2. Tangential Rayleigh model (default)
  D = dampingMassCoef * M + dampingStiffnessCoef * K(q)
     where K is the tangential stiffness matrix IN THE CURRENT CONFIGURATION q

Mode 1. is computationally faster as it does not require LU updates at every 
timestep. However, with reduced geometrically nonlinear FEM simulations, 
the slowdown is very small as the complexity of the timestep is O(r^4) in 
either case, due to the evaluation of reduced internal forces. Mode 1. is 
useful, for example, if you want to timestep linear modal analysis 
simulations (stiffness matrix is constant in that case).

Mode 2. gives a better damping model for large deformations.

In order to use this class, you need to set the timestep very small, or 
else the explicit integrator will go unstable. Roughly speaking, the timestep 
must resolve the highest frequency present in your simulation. With linear 
modal analysis, the spectrum is truncated, so your timestep, while small, 
need not be extremely small. With nonlinear reduced simulations as per [1], 
very high frequencies can be present in the simulation, so you either need to 
take very small explicit timesteps, or use the ImplicitNewmark class.

*/

#ifndef _CENTRALDIFFERENCESDENSE_H_
#define _CENTRALDIFFERENCESDENSE_H_

#include "integratorBaseDense.h"

class CentralDifferencesDense : public virtual IntegratorBaseDense
{
public:
  CentralDifferencesDense(int numDOFs, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, int tangentialDampingMode = 1);

  virtual ~CentralDifferencesDense();

  inline virtual void SetTimestep(double timestep) { this->timestep = timestep; UpdateLU(); }

  // performs one timestep of simulation
  virtual int DoTimestep(); 

  inline void SetTangentialDampingMode(int tangentialDampingMode) { this->tangentialDampingMode = tangentialDampingMode; }

protected:
  double * rhs;
  double * LUFactor;
  double * stiffnessMatrix;

  int tangentialDampingMode;
  int UpdateLU();
};

#endif

