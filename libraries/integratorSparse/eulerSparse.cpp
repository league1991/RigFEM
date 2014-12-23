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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrixIO.h"
#include "performanceCounter.h"
#include "insertRows.h"
#include "eulerSparse.h"

EulerSparse::EulerSparse(int r, double timestep, SparseMatrix * massMatrix_, ForceModel * forceModel_, int symplectic_, int numConstrainedDOFs_, int * constrainedDOFs_, double dampingMassCoef): IntegratorBaseSparse(r, timestep, massMatrix_, forceModel_, numConstrainedDOFs_, constrainedDOFs_, dampingMassCoef, 0.0), symplectic(symplectic_)
{
  #ifdef PARDISO
    printf("Creating Pardiso solver for M.\n");
    int positiveDefiniteSolver = 1;
    pardisoSolver = new PardisoSolver(massMatrix, 1, positiveDefiniteSolver);
    int info = pardisoSolver->ComputeCholeskyDecomposition(massMatrix);
    if (info != 0)
    {
      printf("Error: PARDISO solver returned non-zero exit code %d.\n", info);
      exit(1);
    }
    printf("Solver created.\n");
  #endif

  #ifdef SPOOLES
    printf("Creating SPOOLES solver for M.\n");
    spoolesSolver = new SPOOLESSolver(massMatrix);
    printf("Solver created.\n");
  #endif

  #ifdef PCG
    printf("Creating Jacobi solver for M.\n");
    jacobiPreconditionedCGSolver = new CGSolver(massMatrix);
    printf("Solver created.\n");
  #endif
}

EulerSparse::~EulerSparse()
{
  #ifdef PARDISO
    delete(pardisoSolver);
  #endif

  #ifdef SPOOLES
    delete(spoolesSolver);
  #endif

  #ifdef PCG
    delete(jacobiPreconditionedCGSolver);
  #endif
}

// sets the state based on given q, qvel
// automatically computes acceleration assuming zero external force
int EulerSparse::SetState(double * q_, double * qvel_)
{
  memcpy(q, q_, sizeof(double)*r);

  if (qvel_ != NULL)
    memcpy(qvel, qvel_, sizeof(double)*r);

  return 0;
}

int EulerSparse::DoTimestep()
{
  // v_{n+1} = v_n + h * (F_n / m)
  // x_{n+1} = x_n + h * v_{n+1}

  // store current state
  for(int i=0; i<r; i++)
  {
    q_1[i] = q[i]; 
    qvel_1[i] = qvel[i];
    qaccel_1[i] = qaccel[i];
  }

  PerformanceCounter counterForceAssemblyTime;
  forceModel->GetInternalForce(q, internalForces);
  counterForceAssemblyTime.StopCounter();
  forceAssemblyTime = counterForceAssemblyTime.GetElapsedTime();

  // scale internal forces
  for(int i=0; i<r; i++)
    internalForces[i] *= internalForceScalingFactor;

  // damping
  double * dampingForces = buffer;
  massMatrix->MultiplyVector(qvel, dampingForces);
  for(int i=0; i<r; i++)
    dampingForces[i] *= dampingMassCoef;
  dampingMatrix->MultiplyVectorAdd(qvel, dampingForces);

  //printf("C=\n");
  //dampingMatrix->Print();
  //dampingMatrix->Save("C");

  for(int i=0; i<r; i++)
  {
    // set qresidual = F_n, for a subsequent solve M * qdelta = h * F_n
    qresidual[i] = externalForces[i] - internalForces[i] - dampingForces[i];
  }

  PerformanceCounter counterSystemSolveTime;

  // solve: M * qdelta = qresidual

  memset(qdelta, 0.0, sizeof(double)*r);

  #ifdef PARDISO
    int info = pardisoSolver->SolveLinearSystem(qdelta, qresidual);
    char solverString[16] = "PARDISO";
  #endif

  #ifdef SPOOLES
    int info = spoolesSolver->SolveLinearSystem(qdelta, qresidual);
    char solverString[16] = "SPOOLES";
  #endif

  #ifdef PCG
    int info = jacobiPreconditionedCGSolver->SolveLinearSystemWithJacobiPreconditioner(qdelta, qresidual, 1e-6, 10000);
    if (info > 0)
      info = 0;
    char solverString[16] = "PCG";
  #endif

  if (info != 0)
  {
    printf("Error: %s sparse solver returned non-zero exit status %d.\n", solverString, (int)info);
    return 1;
  }

  counterSystemSolveTime.StopCounter();
  systemSolveTime = counterSystemSolveTime.GetElapsedTime();

  // update state
  if (symplectic)
  {
    for(int i=0; i<r; i++)
    {
      qvel[i] += timestep * qdelta[i];
      q[i] += timestep * qvel[i];
    }	
  }
  else
  {
    for(int i=0; i<r; i++)
    {
      q[i] += timestep * qvel[i];
      qvel[i] += timestep * qdelta[i];
    }
  }

  for(int i=0; i<r; i++)
    qaccel[i] = qdelta[i];

  // constrain fixed DOFs
  for(int i=0; i<numConstrainedDOFs; i++)
    q[constrainedDOFs[i]] = qvel[constrainedDOFs[i]] = qaccel[constrainedDOFs[i]] = 0.0;

  return 0;
}

