/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "reducedForceModel" library , Copyright (C) 2007 CMU, 2009 MIT,       *
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

#include <math.h>
#include "reducedForceModelWithHessian.h"
#include "lapack-headers.h"
#include "matrixMacros.h"

ReducedForceModelWithHessian::ReducedForceModelWithHessian()
{
  hessianBuffer = (double*) malloc (sizeof(double) * r * r * (r+1) / 2);
}

ReducedForceModelWithHessian::~ReducedForceModelWithHessian() 
{ 
  free(hessianBuffer);
}

void ReducedForceModelWithHessian::GetStiffnessMatrixCorrection(double * u, double * du, double * dK) 
{
  GetTangentHessianTensor(u, hessianBuffer);
  ContractTensorWithVector(r, hessianBuffer, du, dK);
}

void ReducedForceModelWithHessian::ContractTensorWithVector(int r, double * Hq, double * q, double * A)
{
  // computes A = Hq : q

  int quadraticSize = r*(r+1)/2;

  // multiply Hq and q
  cblas_dgemv(CblasColMajor, CblasTrans,
       r, quadraticSize,
       1.0,
       Hq, r,
       q, 1,
       0.0,
       A, 1);

  for(int j=r-1; j>=0; j--)
    for(int i=r-1; i>=j; i--)
    {
      int lowerTrianglePos = j * r - (j-1) * j / 2 + (i-j);
      A[ELT(r,i,j)] = A[lowerTrianglePos];
      A[ELT(r,j,i)] = A[lowerTrianglePos];
    }
}

void ReducedForceModelWithHessian::TestHessian(int numTests, double qAmplitude)
{
  double * q = (double*) malloc (sizeof(double) * r);
  double * q1 = (double*) malloc (sizeof(double) * r);
  double * dq = (double*) malloc (sizeof(double) * r);
  double * dqeps = (double*) malloc (sizeof(double) * r);
  double * K = (double*) malloc (sizeof(double) * r * r);
  double * K1 = (double*) malloc (sizeof(double) * r * r);
  double * dK = (double*) malloc (sizeof(double) * r * r);
  for(int i=0; i<numTests; i++)
  {
    for(int j=0; j<r; j++)
    {
      q[j] = 2.0 * (-1.0 + 2.0 * rand() / RAND_MAX); 
      dq[j] = 2.0 * (-1.0 + 2.0 * rand() / RAND_MAX);    
    }

    GetTangentStiffnessMatrix(q, K);
    double eps = 1.0;
    do
    {
      for(int j=0; j<r; j++)
      {
        q1[j] = q[j] + eps * dq[j];
        dqeps[j] = eps * dq[j];
      }

      GetStiffnessMatrixCorrection(q, dqeps, dK);
      GetTangentStiffnessMatrix(q1, K1);

      // K + dK - K1 should be O(eps^2)
      double maxEntry = 0.0;
      for(int j=0; j<r*r; j++)
      {
        double entry = fabs( (K1[j] - K[j] - dK[j]) / eps );
        if (entry > maxEntry)
          maxEntry = entry;
      }

      printf("Hessian trial %d: eps=%G: maxEntry=%G\n", i, eps, maxEntry);
      eps *= 0.5;
    }
    while (eps > 1e-15);
  }
  free(dK);
  free(K1);
  free(K);
  free(dqeps);
  free(dq);
  free(q1);
  free(q);
}

