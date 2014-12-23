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
  A base abstract class to timestep dense dynamics
  (i.e., usually model reduction dynamics).

  It is derived into classes that support dense implicit Newmark,
  and dense central differences.

  See also integratorBase.h.
*/

#ifndef _INTEGRATORBASEDENSE_H_
#define _INTEGRATORBASEDENSE_H_

#include "integratorBase.h"
#include "reducedForceModel.h"

class IPIVC;

// this abstract class is derived into: implicit Newmark, central differences
class IntegratorBaseDense : public IntegratorBase
{
public:
  // r is the dimension of the simulation basis
  // mass matrix is a r x r matrix; you need to pass the identity matrix 
  // when using the method from [1]
  // the damping coefficients are tangential Rayleigh damping coefficients, see [2]
  IntegratorBaseDense(int r, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0);

  virtual ~IntegratorBaseDense();

  // you would rarely need to call this (b/c typically set once and for all in the constructor)
  inline void SetReducedForceModel(ReducedForceModel * reducedForceModel) { this->reducedForceModel = reducedForceModel; }

  void SetMassMatrix(double * massMatrix);
  void SetTangentStiffnessMatrixOffset(double * tangentStiffnessMatrixOffset);
  void ClearTangentStiffnessMatrixOffset();

  // dynamic solver is default (i.e. useStaticSolver=false)
  // with the static solver, all dynamic terms are neglected, and the system only computes the static equilibrium under the currently applied external forces
  void UseStaticSolver(bool useStaticSolver);
  virtual void ResetToRest();

  virtual double GetKineticEnergy();
  virtual double GetTotalMass();

  virtual int SetState(double * q, double * qvel=NULL);

  // returns the execution time of the last r x r linear system solve
  inline virtual double GetForceAssemblyTime() { return forceAssemblyTime; }
  inline virtual double GetSystemSolveTime() { return systemSolveTime; }

  // plastic deformations
  void UsePlasticDeformations(int usePlasticDeformations);
  void SetPlasticThreshold(double plasticThreshold2);
  void ClearPlasticDeformations();

protected:
  
  IntegratorBaseDense(int r, double timestep, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0);

  double * massMatrix; // pointer to the reduced mass matrix
  double * dampingMatrix; 
  double * tangentStiffnessMatrix; 

  #ifdef __APPLE__
    #define INTEGER __CLPK_integer
  #else
    #define INTEGER int
  #endif

  double forceAssemblyTime, systemSolveTime;

  ReducedForceModel * reducedForceModel; 

  int r2; // r * r

  IPIVC * IPIV;

  bool useStaticSolver;
  double * tangentStiffnessMatrixOffset;

  // plastic deformations
  int usePlasticDeformations;
  double plasticThreshold2;
  double * plasticfq;
  double * totalfq;
  void SetTotalForces(double * fq);
  void ProcessPlasticDeformations();
};

#endif

