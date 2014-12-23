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
#include <string.h>
#include "lapack-headers.h"
#include "performanceCounter.h"
#include "implicitNewmarkDense.h"
#include "IPIVC.h"

ImplicitNewmarkDense::ImplicitNewmarkDense(int r, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, solverType solver, double dampingMassCoef, double dampingStiffnessCoef, int maxIterations, double epsilon, double NewmarkBeta, double NewmarkGamma): IntegratorBaseDense(r, timestep, massMatrix, reducedForceModel, dampingMassCoef, dampingStiffnessCoef)
{
  this->solver = solver;

  this->maxIterations = maxIterations; // maxIterations = 1 for semi-implicit
  this->epsilon = epsilon; 

  this->NewmarkBeta = NewmarkBeta;
  this->NewmarkGamma = NewmarkGamma;

  symmetricSolver_lwork = 64 * r;
  symmetricSolver_work = (double*) malloc (sizeof(double) * symmetricSolver_lwork);

  UpdateAlphas();
}

ImplicitNewmarkDense::ImplicitNewmarkDense(int r, double timestep, double dampingMassCoef, double dampingStiffnessCoef, double NewmarkBeta, double NewmarkGamma): IntegratorBaseDense(r, timestep, massMatrix, reducedForceModel, dampingMassCoef, dampingStiffnessCoef)
{
  this->NewmarkBeta = NewmarkBeta;
  this->NewmarkGamma = NewmarkGamma;

  UpdateAlphas();
}


ImplicitNewmarkDense::~ImplicitNewmarkDense()
{
  free(symmetricSolver_work);
}

void ImplicitNewmarkDense::UpdateAlphas()
{
  alpha1 = 1.0 / (NewmarkBeta * timestep * timestep);
  alpha2 = 1.0 / (NewmarkBeta * timestep);
  alpha3 = (1.0 - 2.0 * NewmarkBeta) / (2.0 * NewmarkBeta);
  alpha4 = NewmarkGamma / (NewmarkBeta * timestep);
  alpha5 = 1 - NewmarkGamma/NewmarkBeta;
  alpha6 = (1.0 - NewmarkGamma / (2.0 * NewmarkBeta)) * timestep;
}

int ImplicitNewmarkDense::DoTimestep()
{
  int numIter = 0;

  double error0 = 0; // error after the first step
  double errorQuotient;

  // store current amplitudes and set initial guesses for qaccel, qvel
  // note: these guesses will later be overriden; they are only used to construct the right-hand-side vector (multiplication with M and C)
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

    PerformanceCounter counterForceAssemblyTime;
    reducedForceModel->GetForceAndMatrix(q, internalForces, tangentStiffnessMatrix);
    counterForceAssemblyTime.StopCounter();
    forceAssemblyTime = counterForceAssemblyTime.GetElapsedTime();

    if (plasticfq != NULL)
    {
      SetTotalForces(internalForces);
      for(int i=0; i<r; i++)
        internalForces[i] -= plasticfq[i];
    }

    // scale internal forces
    for(i=0; i<r; i++)
      internalForces[i] *= internalForceScalingFactor;

/*
    printf("internalForceScalingFactor = %G\n", internalForceScalingFactor);
    printf("q:\n");
    for(int i=0; i<r; i++)
      printf("%G ", q[i]);
    printf("\n");

    printf("Internal forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", internalForces[i]);
    printf("\n");
*/

    for(i=0; i<r2; i++)
      tangentStiffnessMatrix[i] *= internalForceScalingFactor;

    for(i=0; i<r2; i++)
      tangentStiffnessMatrix[i] += tangentStiffnessMatrixOffset[i];

/*
    printf("Tangent stiffness matrix:\n");
    for(int i=0; i<r; i++)
    {
      for(int j=0; j<r; j++)
        printf("%.15f ", tangentStiffnessMatrix[r * j + i]);
      printf("\n");
    }
    printf("Tangent stiffness matrix offset:\n");
    for(int i=0; i<r; i++)
    {
      for(int j=0; j<r; j++)
        printf("%.15f ", tangentStiffnessMatrixOffset[r * j + i]);
      printf("\n");
    }
    printf("----\n");
*/

    //WriteMatrixToDisk_("Kr", r, r, tangentStiffnessMatrix);
    //WriteMatrixToDisk_("Mr", r, r, massMatrix);
    //exit(1);

    memset(qresidual, 0, sizeof(double) * r);

    if (useStaticSolver)
    {
      // no operation
    }
    else
    {
      // build effective stiffness: add mass matrix and damping matrix to tangentStiffnessMatrix
      for(i=0; i<r2; i++)
      {
        dampingMatrix[i] = dampingMassCoef * massMatrix[i] + dampingStiffnessCoef * tangentStiffnessMatrix[i];
        tangentStiffnessMatrix[i] += alpha4 * dampingMatrix[i];
        //tangentStiffnessMatrix[i] += alpha3 * massMatrix[i] + gamma * alpha1 * dampingMatrix[i]; // static Rayleigh damping

        // add mass matrix to the effective stiffness matrix
        tangentStiffnessMatrix[i] += alpha1 * massMatrix[i];
      }

      // compute force residual, store it into aux variable qresidual
      // qresidual = M * qaccel + C * qvel - externalForces + internalForces

      // M * qaccel
      cblas_dgemv(CblasColMajor,CblasNoTrans,
        r,r,1.0,massMatrix,r,qaccel,1,0.0,qresidual,1);

      // += C * qvel
      cblas_dgemv(CblasColMajor,CblasNoTrans,
        r,r,1.0,dampingMatrix,r,qvel,1,1.0,qresidual,1);
    }

    // add externalForces, internalForces
    for(i=0; i<r; i++)
    {
      qresidual[i] += internalForces[i] - externalForces[i];
      qresidual[i] *= -1;
      qdelta[i] = qresidual[i];
    }

/*
    printf("internalForceScalingFactor = %G\n", internalForceScalingFactor);

    printf("internal forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", internalForces[i]);
    printf("\n");

    printf("external forces:\n");
    for(int i=0; i<r; i++)
      printf("%G ", externalForces[i]);
    printf("\n");

    printf("mass matrix:\n");
    for(int i=0; i<r*r; i++)
      printf("%G ", massMatrix[i]);
    printf("\n");

    printf("damping matrix:\n");
    for(int i=0; i<r*r; i++)
      printf("%G ", dampingMatrix[i]);
    printf("\n");

    printf("effective stiffness matrix:\n");
    for(int i=0; i<r*r; i++)
      printf("%G ", tangentStiffnessMatrix[i]);
    printf("\n");

    printf("matrix rhs:\n");
    for(int i=0; i<r; i++)
      printf("%G ", qdelta[i]);
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
      // rel error wrt to initial error before performing this iteration
      errorQuotient = error / error0; 
    }

    if ((errorQuotient < epsilon * epsilon) || (error == 0))
    {
      break;
    }

    // solve (effective stiffness) * qdelta = qresidual
    PerformanceCounter counterSystemSolveTime;
    //counterSystemSolveTime.StartCounter(); // it starts automatically in constructor

    switch (solver)
    {
      case generalMatrixSolver:
      {
        INTEGER N = r;
        INTEGER NRHS = 1;
        double * A = tangentStiffnessMatrix;
        INTEGER LDA = r;
        double * B = qdelta;
        INTEGER LDB = r;
        INTEGER INFO;

        #ifdef __APPLE__
          #define DGESV dgesv_
        #else
          #define DGESV dgesv
        #endif

        DGESV ( &N, &NRHS, A, &LDA, IPIV->GetBuf(), B, &LDB, &INFO );

        if (INFO != 0)
        {
          printf("Error: Gaussian elimination solver returned non-zero exit status %d.\n",(int)INFO);
          return 1;
        }
      }
      break;

      case symmetricMatrixSolver:
      {
        // call dsysv ( uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
  
        #ifdef __APPLE__
          #define DSYSV dsysv_
        #else
          #define DSYSV dsysv
        #endif

        char uplo = 'U';
        INTEGER nrhs = 1;
        INTEGER info;
        INTEGER R = r;

        INTEGER symmetricSolver_lworkI = symmetricSolver_lwork;
        DSYSV ( &uplo, &R, &nrhs, tangentStiffnessMatrix, &R, IPIV->GetBuf(), qdelta, &R, symmetricSolver_work, &symmetricSolver_lworkI, &info);

        if (info != 0)
        {
          printf("Error: Symmetric indefinite solver returned non-zero exit status %d.\n",(int)info);
          return 1;
        }
      }
      break;

      case positiveDefiniteMatrixSolver:
      {
        // call dposv ( uplo, n, nrhs, a, lda, b, ldb, info)

        #ifdef __APPLE__
          #define DPOSV dposv_
        #else
          #define DPOSV dposv
        #endif
  
        char uplo = 'U';
        INTEGER nrhs = 1;
        INTEGER info = 0;
        INTEGER R = r;

        DPOSV ( &uplo, &R, &nrhs, tangentStiffnessMatrix, &R, qdelta, &R, &info);

        if (info != 0)
        {
          printf("Error: Positive-definite Cholesky solver returned non-zero exit status %d.\n",(int)info);
          return 1;
        }

      }
      break;

      default:
        printf("Error: reduced integration solver not specified.\n");
        return 1;
      break;
    }
    counterSystemSolveTime.StopCounter();
    systemSolveTime = counterSystemSolveTime.GetElapsedTime();

/*
    printf("qdelta:\n");
    for(int i=0; i<r; i++)
      printf("%G ", qdelta[i]);
    printf("\n");
*/

    // update state
    for(i=0; i<r; i++)
    {
      q[i] += qdelta[i];
      qaccel[i] = alpha1 * (q[i] - q_1[i]) - alpha2 * qvel_1[i] - alpha3 * qaccel_1[i];
      qvel[i] = alpha4 * (q[i] - q_1[i]) + alpha5 * qvel_1[i] + alpha6 * qaccel_1[i];
    }

    numIter++;
  }
  while (numIter < maxIterations);

  ProcessPlasticDeformations();

/*
  printf("Num iterations performed: %d (maxIterations=%d)\n", numIter, maxIterations);
  if ((numIter >= maxIterations) && (maxIterations > 1))
  {
    printf("Warning: method did not converge in max number of iterations.\n");
  }
*/

  return 0;
}

