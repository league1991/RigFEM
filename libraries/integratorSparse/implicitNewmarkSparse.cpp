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
#include "implicitNewmarkSparse.h"

ImplicitNewmarkSparse::ImplicitNewmarkSparse(int r, double timestep, SparseMatrix * massMatrix_, ForceModel * forceModel_, int positiveDefiniteSolver_, int numConstrainedDOFs_, int * constrainedDOFs_, double dampingMassCoef, double dampingStiffnessCoef, int maxIterations, double epsilon, double NewmarkBeta, double NewmarkGamma, int numSolverThreads_): IntegratorBaseSparse(r, timestep, massMatrix_, forceModel_, numConstrainedDOFs_, constrainedDOFs_, dampingMassCoef, dampingStiffnessCoef), positiveDefiniteSolver(positiveDefiniteSolver_), numSolverThreads(numSolverThreads_)
{
  this->maxIterations = maxIterations; // maxIterations = 1 for semi-implicit
  this->epsilon = epsilon; 
  this->NewmarkBeta = NewmarkBeta;
  this->NewmarkGamma = NewmarkGamma;

  useStaticSolver = false;

  UpdateAlphas();

  forceModel->GetTangentStiffnessMatrixTopology(&tangentStiffnessMatrix);

  if (tangentStiffnessMatrix->Getn() != massMatrix->Getn())
  {
    printf("Error: the provided mass matrix does not have correct size. Mass matrix: %d x %d. Stiffness matrix: %d x %d.\n", massMatrix->Getn(), massMatrix->Getn(), tangentStiffnessMatrix->Getn(), tangentStiffnessMatrix->Getn());
    exit(1);
  }

  rayleighDampingMatrix = new SparseMatrix(*tangentStiffnessMatrix);

  rayleighDampingMatrix->BuildSubMatrixIndices(*massMatrix);
  tangentStiffnessMatrix->BuildSubMatrixIndices(*massMatrix);
  tangentStiffnessMatrix->BuildSubMatrixIndices(*dampingMatrix, 1);

  if (tangentStiffnessMatrix->GetNumRows() != massMatrix->GetNumRows())
  {
    printf("Error: mass matrix and stiffness matrix don't have same dimensions.\n");
    exit(1);
  }

  bufferConstrained = (double*) malloc (sizeof(double) * (r - numConstrainedDOFs));

  systemMatrix = new SparseMatrix(*tangentStiffnessMatrix);
  systemMatrix->RemoveRowsColumns(numConstrainedDOFs, constrainedDOFs);
  systemMatrix->BuildSuperMatrixIndices(numConstrainedDOFs, constrainedDOFs, tangentStiffnessMatrix);

  #ifdef PARDISO
    printf("Creating Pardiso solver. Positive-definite solver: %d. Num threads: %d\n", positiveDefiniteSolver, numSolverThreads);
    pardisoSolver = new PardisoSolver(systemMatrix, numSolverThreads, positiveDefiniteSolver);
  #endif

  #ifdef PCG
    jacobiPreconditionedCGSolver = new CGSolver(systemMatrix);
  #endif
}

ImplicitNewmarkSparse::~ImplicitNewmarkSparse()
{
  delete(tangentStiffnessMatrix);
  delete(rayleighDampingMatrix);
  delete(systemMatrix);
  free(bufferConstrained);
  #ifdef PARDISO
    delete(pardisoSolver);
  #endif
}

void ImplicitNewmarkSparse::SetDampingMatrix(SparseMatrix * dampingMatrix)
{
  IntegratorBaseSparse::SetDampingMatrix(dampingMatrix);
  tangentStiffnessMatrix->BuildSubMatrixIndices(*dampingMatrix, 1);
}

void ImplicitNewmarkSparse::UpdateAlphas()
{
  alpha1 = 1.0 / (NewmarkBeta * timestep * timestep);
  alpha2 = 1.0 / (NewmarkBeta * timestep);
  alpha3 = (1.0 - 2.0 * NewmarkBeta) / (2.0 * NewmarkBeta);
  alpha4 = NewmarkGamma / (NewmarkBeta * timestep);
  alpha5 = 1 - NewmarkGamma/NewmarkBeta;
  alpha6 = (1.0 - NewmarkGamma / (2.0 * NewmarkBeta)) * timestep;
}

// sets the state based on given q, qvel
// automatically computes acceleration assuming zero external force
int ImplicitNewmarkSparse::SetState(double * q_, double * qvel_)
{
  memcpy(q, q_, sizeof(double)*r);

  if (qvel_ != NULL)
    memcpy(qvel, qvel_, sizeof(double)*r);
  else
    memset(qvel, 0, sizeof(double)*r);

  for(int i=0; i<numConstrainedDOFs; i++)
    q[constrainedDOFs[i]] = qvel[constrainedDOFs[i]] = 0.0;

  // M * qaccel + C * qvel + R(q) = P_0 
  // R(q) = P_0 = 0
  // i.e. M * qaccel = - C * qvel - R(q)

  forceModel->GetForceAndMatrix(q, internalForces, tangentStiffnessMatrix);

  *rayleighDampingMatrix = dampingStiffnessCoef * (*tangentStiffnessMatrix);
  rayleighDampingMatrix->AddSubMatrix(dampingMassCoef, *massMatrix);

  // buffer = C * qvel
  rayleighDampingMatrix->MultiplyVector(qvel, buffer);
  dampingMatrix->MultiplyVectorAdd(qvel, buffer);

  for(int i=0; i<r; i++)
    buffer[i] = -buffer[i] - internalForces[i];

  // solve M * qaccel = buffer
  RemoveRows(r, bufferConstrained, buffer, numConstrainedDOFs, constrainedDOFs);

  // use tangentStiffnessMatrix as the buffer place
  tangentStiffnessMatrix->ResetToZero();
  tangentStiffnessMatrix->AddSubMatrix(1.0, *massMatrix);
  tangentStiffnessMatrix->AddSubMatrix(1.0, *dampingMatrix, 1);
  systemMatrix->AssignSuperMatrix(tangentStiffnessMatrix); // must go via a matrix with tangentStiffnessMatrix's topology, because the AssignSuperMatrix indices were computed with respect to such topology

  memset(buffer, 0, sizeof(double) * r);

  #ifdef SPOOLES
    SPOOLESSolver solver(systemMatrix);
    int info = solver.SolveLinearSystem(buffer, bufferConstrained);
    char solverString[16] = "SPOOLES";
  #endif

  //massMatrix->Save("M");
  //systemMatrix->Save("A");

  #ifdef PARDISO
    pardisoSolver->ComputeCholeskyDecomposition(systemMatrix);
    int info = pardisoSolver->SolveLinearSystem(buffer, bufferConstrained);
    char solverString[16] = "PARDISO";
  #endif

  #ifdef PCG
    int info = jacobiPreconditionedCGSolver->SolveLinearSystemWithJacobiPreconditioner(buffer, bufferConstrained, 1e-6, 10000);
    if (info > 0)
      info = 0;
    char solverString[16] = "PCG";
  #endif

  if (info != 0)
  {
    printf("Error: %s sparse solver returned non-zero exit status %d.\n", solverString, (int)info);
    return 1;
  }
  
  InsertRows(r, buffer, qaccel, numConstrainedDOFs, constrainedDOFs);

  return 0;
}
 
int ImplicitNewmarkSparse::DoTimestep()
{
  int numIter = 0;

  double error0 = 0; // error after the first step
  double errorQuotient;

  // store current amplitudes and set initial guesses for qaccel, qvel
  for(int i=0; i<r; i++)
  {
    q_1[i] = q[i]; 
    qvel_1[i] = qvel[i];
    qaccel_1[i] = qaccel[i];

    qaccel[i] = alpha1 * (q[i] - q_1[i]) - alpha2 * qvel_1[i] - alpha3 * qaccel_1[i];
    qvel[i] = alpha4 * (q[i] - q_1[i]) + alpha5 * qvel_1[i] + alpha6 * qaccel_1[i];
  }

  do
  {
    int i;

/*
    printf("q:\n");
    for(int i=0; i<r; i++)
      printf("%G ", q[i]);
    printf("\n");

    printf("Internal forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", internalForces[i]);
    printf("\n");
*/

    PerformanceCounter counterForceAssemblyTime;
    forceModel->GetForceAndMatrix(q, internalForces, tangentStiffnessMatrix);
    counterForceAssemblyTime.StopCounter();
    forceAssemblyTime = counterForceAssemblyTime.GetElapsedTime();

    //tangentStiffnessMatrix->Print();
    //tangentStiffnessMatrix->Save("K");

    // scale internal forces
    for(i=0; i<r; i++)
      internalForces[i] *= internalForceScalingFactor;

    *tangentStiffnessMatrix *= internalForceScalingFactor;

    memset(qresidual, 0, sizeof(double) * r);

    if (useStaticSolver)
    {
      // no operation
    }
    else
    {
      // build effective stiffness: add mass matrix and damping matrix to tangentStiffnessMatrix
      tangentStiffnessMatrix->ScalarMultiply(dampingStiffnessCoef, rayleighDampingMatrix);
      rayleighDampingMatrix->AddSubMatrix(dampingMassCoef, *massMatrix);

      rayleighDampingMatrix->ScalarMultiplyAdd(alpha4, tangentStiffnessMatrix);
      //*tangentStiffnessMatrix += alpha4 * *rayleighDampingMatrix;
      tangentStiffnessMatrix->AddSubMatrix(alpha4, *dampingMatrix, 1);

      tangentStiffnessMatrix->AddSubMatrix(alpha1, *massMatrix);
      
      // compute force residual, store it into aux variable qresidual
      // qresidual = M * qaccel + C * qvel - externalForces + internalForces

      massMatrix->MultiplyVector(qaccel, qresidual);
      rayleighDampingMatrix->MultiplyVectorAdd(qvel, qresidual);
      dampingMatrix->MultiplyVectorAdd(qvel, qresidual);
    }

    // add externalForces, internalForces
    for(i=0; i<r; i++)
    {
      qresidual[i] += internalForces[i] - externalForces[i];
      qresidual[i] *= -1;
      qdelta[i] = qresidual[i];
    }

/*
    printf("internal forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", internalForces[i]);
    printf("\n");

    printf("external forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", externalForces[i]);
    printf("\n");

    printf("residual:\n");
    for(int i=0; i<r; i++)
      printf("%G ", -qresidual[i]);
    printf("\n");
*/

    double error = 0;
    for(i=0; i<r; i++)
      error += qresidual[i] * qresidual[i];

    // on the first iteration, compute initial error
    if (numIter == 0) 
    {
      error0 = error;
      errorQuotient = 1.0;
    }
    else
    {
      // error divided by the initial error, before performing this iteration
      errorQuotient = error / error0; 
    }

    if (errorQuotient < epsilon * epsilon)
    {
      break;
    }

    //tangentStiffnessMatrix->Save("Keff");
    RemoveRows(r, bufferConstrained, qdelta, numConstrainedDOFs, constrainedDOFs);
    systemMatrix->AssignSuperMatrix(tangentStiffnessMatrix);

    // solve: systemMatrix * buffer = bufferConstrained

    PerformanceCounter counterSystemSolveTime;
    memset(buffer, 0, sizeof(double) * r);

    #ifdef SPOOLES
      SPOOLESSolver solver(systemMatrix);
      int info = solver.SolveLinearSystem(buffer, bufferConstrained);
      char solverString[16] = "SPOOLES";
    #endif

    #ifdef PARDISO
      int info = pardisoSolver->ComputeCholeskyDecomposition(systemMatrix);
      if (info == 0)
        info = pardisoSolver->SolveLinearSystem(buffer, bufferConstrained);
      char solverString[16] = "PARDISO";
    #endif

    #ifdef PCG
      int info = jacobiPreconditionedCGSolver->SolveLinearSystemWithJacobiPreconditioner(buffer, bufferConstrained, 1e-6, 10000);
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

    InsertRows(r, buffer, qdelta, numConstrainedDOFs, constrainedDOFs);

/*
    printf("qdelta:\n");
    for(int i=0; i<r; i++)
      printf("%G ", qdelta[i]);
    printf("\n");
    exit(1);
*/
    // update state
    for(i=0; i<r; i++)
    {
      q[i] += qdelta[i];
      qaccel[i] = alpha1 * (q[i] - q_1[i]) - alpha2 * qvel_1[i] - alpha3 * qaccel_1[i];
      qvel[i] = alpha4 * (q[i] - q_1[i]) + alpha5 * qvel_1[i] + alpha6 * qaccel_1[i];
    }

    for(int i=0; i<numConstrainedDOFs; i++)
      q[constrainedDOFs[i]] = qvel[constrainedDOFs[i]] = qaccel[constrainedDOFs[i]] = 0.0;

    numIter++;
  }
  while (numIter < maxIterations);

/*
  printf("qvel:\n");
  for(int i=0; i<r; i++)
    printf("%G ", qvel[i]);
  printf("\n");

  printf("qaccel:\n");
  for(int i=0; i<r; i++)
    printf("%G ", qaccel[i]);
  printf("\n");
*/

  //printf("Num iterations performed: %d\n",numIter);
  //if ((numIter >= maxIterations) && (maxIterations > 1))
  //{
    //printf("Warning: method did not converge in max number of iterations.\n");
  //}

  return 0;
}

void ImplicitNewmarkSparse::UseStaticSolver(bool useStaticSolver_)
{ 
  useStaticSolver = useStaticSolver_;

  if (!useStaticSolver) 
  {
    memset(qvel, 0, sizeof(double) * r);
    memset(qaccel, 0, sizeof(double) * r);
    memset(qvel_1, 0, sizeof(double) * r); 
    memset(qaccel_1, 0, sizeof(double) * r);
    memcpy(q_1, q, sizeof(double) * r);
  }
} 

