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

#include "reducedLinearForceModel.h"
#include "lapack-headers.h"

ReducedLinearForceModel::ReducedLinearForceModel(int r, double * stiffnessMatrix) : ReducedForceModel()
{
  this->r = r;
  this->stiffnessMatrix = (double*) malloc (sizeof(double)*r*r);
  memcpy(this->stiffnessMatrix, stiffnessMatrix, sizeof(double)*r*r);
}

ReducedLinearForceModel::~ReducedLinearForceModel()
{
  free(stiffnessMatrix);
}

void ReducedLinearForceModel::GetInternalForce(double * q, double * internalForces)
{
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasNoTrans;
  int m = r;
  int n = r;
  double alpha = 1;
  double * a = stiffnessMatrix;
  int lda = r;
  double * x = q;
  int incx = 1;
  double beta = 0;
  double * y = internalForces; 
  int incy = 1;

  cblas_dgemv(order, trans, m, n, alpha, a, lda, x, incx,
                  beta, y, incy);
}

void ReducedLinearForceModel::GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix)
{
  memcpy(tangentStiffnessMatrix,stiffnessMatrix,sizeof(double)*r*r);
}

