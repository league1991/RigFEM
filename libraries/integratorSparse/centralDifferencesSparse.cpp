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
#include "performanceCounter.h"
#include "insertRows.h"
#include "centralDifferencesSparse.h"

CentralDifferencesSparse::CentralDifferencesSparse(int numDOFs, double timestep, SparseMatrix * massMatrix_, ForceModel * forceModel_, int numConstrainedDOFs, int * constrainedDOFs, double dampingMassCoef, double dampingStiffnessCoef, int tangentialDampingMode_, int numSolverThreads_): IntegratorBaseSparse(numDOFs, timestep, massMatrix_, forceModel_, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef), tangentialDampingMode(tangentialDampingMode_), numSolverThreads(numSolverThreads_), timestepIndex(0)
{
  rhs = (double*) malloc (sizeof(double) * r);
  rhsConstrained = (double*) malloc (sizeof(double) * (r - numConstrainedDOFs));

  forceModel->GetTangentStiffnessMatrixTopology(&tangentStiffnessMatrix);
  rayleighDampingMatrix = new SparseMatrix(*tangentStiffnessMatrix);
  rayleighDampingMatrix->BuildSubMatrixIndices(*massMatrix);
  tangentStiffnessMatrix->BuildSubMatrixIndices(*massMatrix);

  if (tangentStiffnessMatrix->GetNumRows() != massMatrix->GetNumRows())
  {
    printf("Error: mass matrix and stiffness matrix don't have same dimensions.\n");
    exit(1);
  }

  systemMatrix = new SparseMatrix(*tangentStiffnessMatrix);
  systemMatrix->RemoveRowsColumns(numConstrainedDOFs, constrainedDOFs);
  systemMatrix->BuildSuperMatrixIndices(numConstrainedDOFs, constrainedDOFs, tangentStiffnessMatrix);

  #ifdef PARDISO
    printf("Creating Pardiso solver for central differences.\n");
    int positiveDefiniteSolver = 0;
    pardisoSolver = new PardisoSolver(systemMatrix, numSolverThreads, positiveDefiniteSolver);
  #endif

  #ifdef SPOOLES
    spoolesSolver = NULL;
  #endif

  #ifdef PCG
    printf("Creating Jacobi solver for central differences.\n");
    jacobiPreconditionedCGSolver = new CGSolver(systemMatrix);
  #endif

  DecomposeSystemMatrix();
}

CentralDifferencesSparse::~CentralDifferencesSparse()
{
  #ifdef PARDISO
    delete(pardisoSolver);
  #endif

  #ifdef SPOOLES
    if (spoolesSolver != NULL)
      delete(spoolesSolver);
  #endif

  #ifdef PCG
    delete(jacobiPreconditionedCGSolver);
  #endif

  delete(systemMatrix);
  delete(tangentStiffnessMatrix);
  delete(rayleighDampingMatrix);
  free(rhs);
  free(rhsConstrained);
}

void CentralDifferencesSparse::DecomposeSystemMatrix()
{
  //printf("*** Central differences: decomposing the system matrix.\n");
  // construct damping matrix
  // rayleigh damping matrix = dampingMasscoef * massMatrix + dampingStiffnessCoef * stiffness matrix
  forceModel->GetTangentStiffnessMatrix(q, tangentStiffnessMatrix);
  tangentStiffnessMatrix->ScalarMultiply(internalForceScalingFactor);

  tangentStiffnessMatrix->ScalarMultiply(dampingStiffnessCoef, rayleighDampingMatrix);
  rayleighDampingMatrix->AddSubMatrix(dampingMassCoef, *massMatrix);

  // system matrix = mass matrix + 0.5 * timestep * damping matrix (and remove constrained rows and columns)
  rayleighDampingMatrix->ScalarMultiply(0.5 * timestep, tangentStiffnessMatrix);
  tangentStiffnessMatrix->AddSubMatrix(1.0, *massMatrix);
  systemMatrix->AssignSuperMatrix(tangentStiffnessMatrix);

  //systemMatrix->SaveToMatlabFormat("system.mat");
  
  #ifdef PARDISO
    int info = pardisoSolver->ComputeCholeskyDecomposition(systemMatrix);
    if (info != 0)
    {
      printf("Error: PARDISO solver returned non-zero exit code %d.\n", info);
      exit(1);
    }
  #endif

  #ifdef SPOOLES
    printf("Creating SPOOLES solver for central differences.\n");
    if (spoolesSolver != NULL)
      delete(spoolesSolver);
    if (numSolverThreads > 1)
      spoolesSolver = new SPOOLESSolverMT(systemMatrix, numSolverThreads);
    else
      spoolesSolver = new SPOOLESSolver(systemMatrix);
  #endif
}

int CentralDifferencesSparse::DoTimestep()
{
  PerformanceCounter counterForceAssemblyTime;
    forceModel->GetInternalForce(q, internalForces);
    for (int i=0; i<r; i++)
      internalForces[i] *= internalForceScalingFactor;
  counterForceAssemblyTime.StopCounter();
  forceAssemblyTime = counterForceAssemblyTime.GetElapsedTime();

  if (tangentialDampingMode > 0)
    if (timestepIndex % tangentialDampingMode == 0)
      DecomposeSystemMatrix(); // this routines also updates the damping and system matrices
  
  // update equation is (see WRIGGERS P.: Computational Contact Mechanics. John Wiley & Sons, Ltd., 2002., page 275) :
  //
  // (M + dt / 2 * C) * q(t+1) = (dt)^2 * (fext(t) - fint(q(t))) + dt / 2 * C * q(t-1) + M * (2q(t) - q(t-1))
  //
  // (M + dt / 2 * C) * (q(t+1) - q(t)) = (dt)^2 * (fext(t) - fint(q(t))) + dt / 2 * C * (q(t-1) - q(t)) + M * (q(t) - q(t-1)) 

  // fext are the external forces
  // fint is the vector of internal forces

  // compute rhs = (dt)^2 * (fext - fint(q(t))) + dt / 2 * C * (q(t-1) - q(t)) + M * (q(t) - q(t-1))
  // first, compute rhs = M * (q - q_1)
  for (int i=0; i<r; i++)
    buffer[i] = q[i] - q_1[i];
  massMatrix->MultiplyVector(buffer, rhs);
  
  // rhs += dt / 2 * dampingMatrix * (q_{n-1} - q_n)
  for (int i=0; i<r; i++)
    qdelta[i] = q_1[i] - q[i];
  rayleighDampingMatrix->MultiplyVector(qdelta, buffer);
  for (int i=0; i<r; i++)
    rhs[i] += 0.5 * timestep * buffer[i];

  // rhs += dt * dt * (fext - fint(q(t))) 
  double timestep2 = timestep * timestep;
  for (int i=0; i<r; i++)
    rhs[i] += timestep2 * (externalForces[i] - internalForces[i]);

  // now rhs contains the correct value

  RemoveRows(r, rhsConstrained, rhs, numConstrainedDOFs, constrainedDOFs);

  PerformanceCounter counterSystemSolveTime;

  memset(buffer, 0, sizeof(double) * r);

  #ifdef SPOOLES
    int info = spoolesSolver->SolveLinearSystem(buffer, rhsConstrained);
    char solverString[16] = "SPOOLES";
  #endif

  #ifdef PARDISO
    int info = pardisoSolver->SolveLinearSystem(buffer, rhsConstrained);
    char solverString[16] = "PARDISO";
  #endif
  
  #ifdef PCG
    int info = jacobiPreconditionedCGSolver->SolveLinearSystemWithJacobiPreconditioner(buffer, rhsConstrained, 1e-6, 10000);
    if (info > 0)
      info = 0;
    char solverString[16] = "PCG";
  #endif

  InsertRows(r, buffer, qdelta, numConstrainedDOFs, constrainedDOFs);

  counterSystemSolveTime.StopCounter();
  systemSolveTime = counterSystemSolveTime.GetElapsedTime();

  if (info != 0)
  {
    printf("Error: %s sparse solver returned non-zero exit status %d.\n", solverString, (int)info);
    return 1;
  }

  // the new value of q is now in buffer
  // update velocity, and previous and current positions
  for (int i=0; i<r; i++)
  {
    q_1[i] = q[i];
    qvel[i] = qdelta[i] / timestep;
    qaccel[i] = (qvel[i] - qvel_1[i]) / timestep;
    qvel_1[i] = qvel[i];
    qaccel_1[i] = qaccel[i];
    q[i] += qdelta[i];
  }

  timestepIndex++;

  return 0;
}

// sets the state based on given q, qvel
// automatically computes acceleration assuming zero external force
int CentralDifferencesSparse::SetState(double * q_, double * qvel_)
{
  memcpy(q, q_, sizeof(double) * r);
  memcpy(q_1, q_, sizeof(double) * r);

  if (qvel_ != NULL)
  {
    memcpy(qvel, qvel_, sizeof(double) * r);
    memcpy(qvel_1, qvel_, sizeof(double) * r);
    for(int i=0; i<r; i++)
      q_1[i] = q[i] - timestep * qvel[i];
  }

  return 0;
}

void CentralDifferencesSparse::SetInternalForceScalingFactor(double internalForceScalingFactor)
{
  IntegratorBaseSparse::SetInternalForceScalingFactor(internalForceScalingFactor);
  DecomposeSystemMatrix();
}

void CentralDifferencesSparse::ResetToRest()
{
  IntegratorBaseSparse::ResetToRest();
  timestepIndex = 0;
}

