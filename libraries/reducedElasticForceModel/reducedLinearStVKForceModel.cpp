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
#include "reducedLinearStVKForceModel.h"

ReducedLinearStVKForceModel::ReducedLinearStVKForceModel(StVKReducedStiffnessMatrix * stVKReducedStiffnessMatrix)
{
  r = stVKReducedStiffnessMatrix->Getr();
  K = (double*) malloc (sizeof(double) * r * r);

  double * zero = (double*) calloc (r, sizeof(double));
  stVKReducedStiffnessMatrix->Evaluate(zero,K);
  free(zero);
}

ReducedLinearStVKForceModel::ReducedLinearStVKForceModel(int r, double * K)
{
  this->r = r;
  this->K = (double*) malloc (sizeof(double) * r * r);
  memcpy(this->K, K, sizeof(double) * r * r);
}

void ReducedLinearStVKForceModel::GetInternalForce(double * q, double * internalForces)
{
  // internalForces = K * q
  memset(internalForces, 0, sizeof(double) * r);
  for(int i=0; i<r; i++)
    for(int j=0; j<r; j++)
      internalForces[i] += K[ELT(r,i,j)] * q[j];
}

void ReducedLinearStVKForceModel::GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix)
{
  memcpy(tangentStiffnessMatrix, K, sizeof(double) * r * r);
}

void ReducedLinearStVKForceModel::GetTangentHessianTensor(double * q, double * tangentHessianTensor)
{
  memset(tangentHessianTensor, 0, sizeof(double) * r * r*(r+1)/2);
}

