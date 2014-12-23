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

#include "lapack-headers.h"
#include "integratorBaseDense.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "IPIVC.h"

IntegratorBaseDense::IntegratorBaseDense(int r, double timestep, double * massMatrix, ReducedForceModel * reducedForceModel, double dampingMassCoef, double dampingStiffnessCoef) : IntegratorBase(r, timestep, dampingMassCoef, dampingStiffnessCoef), useStaticSolver(0), usePlasticDeformations(0), plasticThreshold2(DBL_MAX), plasticfq(NULL), totalfq(NULL)
{
  r2 = r * r;
  this->massMatrix = (double*) malloc (sizeof(double) * r2);
  dampingMatrix = (double*) malloc (sizeof(double) * r2);
  tangentStiffnessMatrix = (double*) malloc (sizeof(double) * r2);
  memcpy(this->massMatrix, massMatrix,sizeof(double) * r2);

  forceAssemblyTime = systemSolveTime = 0.0;

  this->reducedForceModel = reducedForceModel;

  tangentStiffnessMatrixOffset = (double*) calloc (r2, sizeof(double));
  IPIV = new IPIVC(r);

  ResetToRest();

  memset(externalForces,0,sizeof(double) * r);
  memset(internalForces,0,sizeof(double) * r);
}

IntegratorBaseDense::IntegratorBaseDense(int r, double timestep, double dampingMassCoef, double dampingStiffnessCoef): IntegratorBase(r, timestep, dampingMassCoef, dampingStiffnessCoef), useStaticSolver(0), usePlasticDeformations(0), plasticThreshold2(DBL_MAX), plasticfq(NULL), totalfq(NULL)
{
  this->massMatrix = NULL;
  dampingMatrix = NULL;
  tangentStiffnessMatrix = NULL;

  forceAssemblyTime = systemSolveTime = 0.0;

  this->reducedForceModel = NULL;

  tangentStiffnessMatrixOffset = NULL;
  IPIV = NULL;

  memset(externalForces,0,sizeof(double) * r);
  memset(internalForces,0,sizeof(double) * r);
}

IntegratorBaseDense::~IntegratorBaseDense()
{
  free(massMatrix);
  free(dampingMatrix);
  free(tangentStiffnessMatrix);
  delete(IPIV);
  free(tangentStiffnessMatrixOffset);
  free(plasticfq);
  free(totalfq);
}

void IntegratorBaseDense::SetMassMatrix(double * massMatrix)
{
  memcpy(this->massMatrix, massMatrix, sizeof(double) * r * r);
}

void IntegratorBaseDense::SetTangentStiffnessMatrixOffset(double * tangentStiffnessMatrixOffset)
{
  memcpy(this->tangentStiffnessMatrixOffset, tangentStiffnessMatrixOffset, sizeof(double) * r * r);
}

void IntegratorBaseDense::ClearTangentStiffnessMatrixOffset()
{
  memset(this->tangentStiffnessMatrixOffset, 0, sizeof(double) * r * r);
}

double IntegratorBaseDense::GetKineticEnergy()
{
  // Wkin = 0.5 * <massMatrix * vel, vel>

  // massMatrix * qvel
  cblas_dgemv(CblasColMajor,CblasNoTrans,r,r,1.0,massMatrix,r,qvel,1,0.0,buffer,1);

  return 0.5 * cblas_ddot(r,qvel,1,buffer,1);
}

double IntegratorBaseDense::GetTotalMass()
{
  int r2 = r*r;
  double totalMass = 0.0;

  for(int i=0; i<r2; i++)
    totalMass += massMatrix[i];

  return totalMass;
}

void IntegratorBaseDense::ResetToRest()
{
  IntegratorBase::ResetToRest();

  if (reducedForceModel != NULL)
    reducedForceModel->ResetToZero();

  ClearPlasticDeformations();
}

int IntegratorBaseDense::SetState(double * q_, double * qvel_)
{
  memcpy(q, q_, sizeof(double)*r);

  if (qvel_ != NULL)
    memcpy(qvel, qvel_, sizeof(double)*r);
  else
    memset(qvel, 0, sizeof(double)*r);

  // M * qaccel + C * qvel + R(q) = P_0 
  // assume P_0 = 0
  // i.e. M * qaccel = - C * qvel - R(q)

  reducedForceModel->GetForceAndMatrix(q, internalForces, tangentStiffnessMatrix);

  int r2 = r * r;
  for(int i=0; i<r2; i++)
    dampingMatrix[i] = dampingMassCoef * massMatrix[i] + dampingStiffnessCoef * tangentStiffnessMatrix[i];

  // buffer = C * qvel
  cblas_dgemv(CblasColMajor, CblasNoTrans,
    r, r, 1.0, dampingMatrix, r, qvel, 1, 0.0, qaccel, 1);

  for(int i=0; i<r; i++)
    qaccel[i] = -qaccel[i] - internalForces[i];

  // solve M * x = qaccel
  // call dposv ( uplo, n, nrhs, a, lda, b, ldb, info)

  #ifdef __APPLE__
    #define DPOSV dposv_  
  #else
    #define DPOSV dposv
  #endif

  memcpy(tangentStiffnessMatrix, massMatrix, sizeof(double) * r2); // must copy mass matrix to another buffer since DPOSV overwrites the system matrix
  char uplo = 'U';
  INTEGER nrhs = 1;
  INTEGER info = 0;
  INTEGER R = r;  
  DPOSV ( &uplo, &R, &nrhs, tangentStiffnessMatrix, &R, qaccel, &R, &info);

  if (info != 0)
  {
    printf("Error: Positive-definite Cholesky solver returned non-zero exit status %d.\n",(int)info);
    return 1;
  }

  return 0;
}

void IntegratorBaseDense::UseStaticSolver(bool useStaticSolver_)
{
  useStaticSolver = useStaticSolver_;

  if (!useStaticSolver)
  {
    memset(qvel,0,sizeof(double)*r);
    memset(qaccel,0,sizeof(double)*r);
    memset(qvel_1,0,sizeof(double)*r);
    memset(qaccel_1,0,sizeof(double)*r);
    memcpy(q_1, q, sizeof(double)*r);
  }
}

void IntegratorBaseDense::UsePlasticDeformations(int usePlasticDeformations_)
{
  usePlasticDeformations = usePlasticDeformations_;
  if (usePlasticDeformations)
  {
    if (plasticfq == NULL)
     plasticfq = (double*) calloc (r, sizeof(double));

    if (totalfq == NULL)
      totalfq = (double*) calloc (r, sizeof(double));
  }
}

void IntegratorBaseDense::SetPlasticThreshold(double plasticThreshold2_)
{
  plasticThreshold2 = plasticThreshold2_;
}

void IntegratorBaseDense::ClearPlasticDeformations()
{
  if (!usePlasticDeformations)
    return;

  for(int i=0; i<r; i++)
  {
    plasticfq[i] = 0.0;
    totalfq[i] = 0.0;
  }
}

void IntegratorBaseDense::ProcessPlasticDeformations()
{
  if (!usePlasticDeformations)
    return;

  double norm2 = 0.0;
  for(int i=0; i<r; i++)
    norm2 += (totalfq[i] - plasticfq[i]) * (totalfq[i] - plasticfq[i]);

  if (norm2 > plasticThreshold2)
  {
    // plastic threshold has been exceeded
    //printf("Plastic threshold exceeded!\n");
    double norm = sqrt(norm2);
    double offset = norm - sqrt(plasticThreshold2);
    for(int i=0; i<r; i++)
      plasticfq[i] += offset / norm * (totalfq[i] - plasticfq[i]);
  }
}

void IntegratorBaseDense::SetTotalForces(double * totalfq_)
{
  memcpy(totalfq, totalfq_, sizeof(double) * r);
}

