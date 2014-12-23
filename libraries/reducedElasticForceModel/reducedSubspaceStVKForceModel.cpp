/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "elasticForceModel" library , Copyright (C) 2007 CMU, 2009 MIT,       *
 *                                                       2014 USC        *
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

#include "matrixMacros.h"
#include "reducedSubspaceStVKForceModel.h"

ReducedSubspaceStVKForceModel::ReducedSubspaceStVKForceModel(
                          StVKInternalForces * stVKInternalForces, 
                          StVKStiffnessMatrix * stVKStiffnessMatrix,
                          ModalMatrix * modalMatrix)
{
  this->stVKInternalForces = stVKInternalForces;
  this->stVKStiffnessMatrix = stVKStiffnessMatrix;
  this->modalMatrix = modalMatrix;
  n = modalMatrix->Getn();
  r = modalMatrix->Getr();
  U = modalMatrix->GetMatrix();

  // allocate buffers
  u = (double*) malloc (sizeof(double) * 3 * n);
  bufferVector = (double*) malloc (sizeof(double) * 3 * n);
  bufferMatrix = (double*) malloc (sizeof(double) * 3 * n * r);

  // generate sparse matrix
  stVKStiffnessMatrix->GetStiffnessMatrixTopology(&sparseMatrix);
}

ReducedSubspaceStVKForceModel::~ReducedSubspaceStVKForceModel()
{
  free(bufferVector);
  free(bufferMatrix);
  free(u);
  delete(sparseMatrix);
}

void ReducedSubspaceStVKForceModel::GetInternalForce(double * q, double * internalForces)
{
  // construct u = U * q
  modalMatrix->AssembleVector(q, u);

  // evaluate forces
  stVKInternalForces->ComputeForces(u, bufferVector);

  // project forces
  modalMatrix->ProjectVector(bufferVector, internalForces);
}

void ReducedSubspaceStVKForceModel::GetTangentStiffnessMatrix(
  double * q, double * tangentStiffnessMatrix)
{
  // construct u = U * q
  modalMatrix->AssembleVector(q, u);
  GetTangentStiffnessMatrixHelper(tangentStiffnessMatrix);
}

void ReducedSubspaceStVKForceModel::GetForceAndMatrix(double * q, 
    double * internalForces, double * tangentStiffnessMatrix)
{
  GetInternalForce(q, internalForces);
  GetTangentStiffnessMatrixHelper(tangentStiffnessMatrix);
}

void ReducedSubspaceStVKForceModel::GetTangentStiffnessMatrixHelper(
  double * tangentStiffnessMatrix)
{
  // evaluate stiffness matrix
  stVKStiffnessMatrix->ComputeStiffnessMatrix(u, sparseMatrix);

  // project matrix
  for(int i=0; i<r; i++)
    sparseMatrix->MultiplyVector(&U[ELT(3*n,0,i)], &bufferMatrix[ELT(3*n,0,i)]);

  modalMatrix->ProjectMatrix(r, bufferMatrix, tangentStiffnessMatrix);
}

