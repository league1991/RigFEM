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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "lapack-headers.h"
#include "matrixMacros.h"
#include "performanceCounter.h"
#include "centralDifferencesDense.h"
#include "IPIVC.h"

#ifdef __APPLE__
  #define DGETRF dgetrf_
  #define DGETRS dgetrs_
#else
  #define DGETRF dgetrf
  #define DGETRS dgetrs
#endif

CentralDifferencesDense::CentralDifferencesDense(int numDOFs, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, double dampingMassCoef, double dampingStiffnessCoef, int tangentialDampingMode_): IntegratorBaseDense(numDOFs, timestep, massMatrix, reducedForceModel, dampingMassCoef, dampingStiffnessCoef), tangentialDampingMode(tangentialDampingMode_)
{
  rhs = (double*) malloc (sizeof(double) * r);
  LUFactor = (double*) malloc (sizeof(double) * r2);
  stiffnessMatrix = (double*) malloc (sizeof(double) * r2);

  UpdateLU();
}

CentralDifferencesDense::~CentralDifferencesDense()
{
  free(stiffnessMatrix);
  free(LUFactor);
  free(rhs);
}

int CentralDifferencesDense::UpdateLU()
{
  if (r == 0)
    return 0;

  reducedForceModel->GetTangentStiffnessMatrix(q, stiffnessMatrix);

  // construct damping matrix
  for(int i=0; i<r2; i++)
    dampingMatrix[i] = dampingMassCoef * massMatrix[i] + dampingStiffnessCoef * stiffnessMatrix[i];

  // do LU decomposition
  for(int i=0; i<r2; i++)
    LUFactor[i] = massMatrix[i] + 0.5 * timestep * dampingMatrix[i];

  INTEGER M = r;
  INTEGER N = r;
  double * A = LUFactor;
  INTEGER LDA = r;
  INTEGER INFO;

  DGETRF(&M, &N, A, &LDA, IPIV->GetBuf(), &INFO);

  if (INFO != 0)
  {
    printf("Warning: LAPACK routine DGETRF returned non-zero exit status %d.\n",(int)INFO);
    return INFO;
  }  

  return 0;
}

int CentralDifferencesDense::DoTimestep()
{
  if (r == 0)
    return 0;

  // the reduced force interpolation
  PerformanceCounter counterForceAssemblyTime;
  reducedForceModel->GetInternalForce(q,internalForces);
  counterForceAssemblyTime.StopCounter();
  forceAssemblyTime = counterForceAssemblyTime.GetElapsedTime();

  if (plasticfq != NULL)
  {
    SetTotalForces(internalForces);
    for(int i=0; i<r; i++)
      internalForces[i] -= plasticfq[i];
  }

  PerformanceCounter counterSystemSolveTime;

  for (int i=0; i<r; i++)
    internalForces[i] *= internalForceScalingFactor;

  if (tangentialDampingMode)
    UpdateLU();
  
  // update equation is:
  //
  // (massMatrix + dt / 2 * dampingMatrix) * q(t+1) = (dt)^2 * (fr - Rr(q(t))) + dt/2 * dampingMatrix * q(t-1) + massMatrix * (2q(t) - q(t-1))

  // LU decomposition of massMatrix + dt / 2 * Dr is available in L,U
  // fr = U^T * f are the reduced external forces
  // Rr is the vector of reduced internal forces

  // update equation follows from Newton's law
  // Mu'' = -Cu' - F_int + F_ext
  // Mu'' + Cu' + R(u) =  F_ext
  // R(u) = F_int

  // here, F_int is the external loading force necessary to sustain a certain deformation
  // it is opposite to the internal forces acting on the body in a given deformation state

  // compute rhs = (dt)^2 * (fr - Rr(q(t))) + dt/2 * dampingMatrix * q(t-1) + massMatrix * (2q(t) - q(t-1))
  // first, compute rhs = massMatrix * (2*q - q_1)
  for (int i=0; i<r; i++)
  {
    rhs[i] = 0;
    for (int j=0; j<r; j++)
      rhs[i] += massMatrix[ELT(r,i,j)] * (2 * q[j] - q_1[j]);
  }

  // rhs += dt / 2 * dampingMatrix * q_{n-1}
  for (int i=0; i<r; i++)
    for (int j=0; j<r; j++)
      rhs[i] += timestep / 2 * dampingMatrix[ELT(r,i,j)] * q_1[j];

  // rhs += dt * dt * (fr - Rr(q(t))) 
  for (int i=0; i<r; i++)
    rhs[i] += timestep * timestep * (externalForces[i] - internalForces[i]);

  // now rhs contains the correct values

  // solve (M~ + dt/2 D~) * qnew = rhs
  // use data from the previously computed LU decomposition
  char trans='N';
  INTEGER nrhs = 1;
  INTEGER INFO;
  INTEGER R = r;
  DGETRS (&trans,&R,&nrhs,LUFactor,&R,IPIV->GetBuf(),rhs,&R,&INFO);

  if (INFO != 0)
  {
    printf("Error: DGETRS returned a non-zero exit code %d.\n", (int)INFO);
    return INFO;
  }

  // the solution qnew is now in rhs
  // update velocity
  // and update previous and current positions
  for (int i=0; i<r; i++)
  {
    qvel[i] = (rhs[i] - q[i]) / timestep;
    q_1[i] = q[i];
    q[i] = rhs[i];
  }

  ProcessPlasticDeformations();

  counterSystemSolveTime.StopCounter();
  systemSolveTime = counterSystemSolveTime.GetElapsedTime();

  return 0;
}

