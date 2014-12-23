/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "forceModelBase" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC *
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
#include <math.h>
#include "forceModel.h"

ForceModel::~ForceModel()
{
}

void ForceModel::GetForceAndMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
  GetInternalForce(u, internalForces);
  GetTangentStiffnessMatrix(u, tangentStiffnessMatrix);
}

void ForceModel::TestStiffnessMatrix(double * q, double * dq)
{
  double * q1 = (double*) malloc (sizeof(double) * r);
  double * dqeps = (double*) malloc (sizeof(double) * r);
  double * internalForce0 = (double*) malloc (sizeof(double) * r);
  double * internalForce1 = (double*) malloc (sizeof(double) * r);
  double * testV = (double*) malloc (sizeof(double) * r);

  SparseMatrix * K0;
  SparseMatrix * K1;
  GetTangentStiffnessMatrixTopology(&K0);
  GetTangentStiffnessMatrixTopology(&K1);
  printf("r: %d K: %d %d\n", r, K0->GetNumRows(), K0->GetNumColumns());

  GetForceAndMatrix(q, internalForce0, K0);
  double KNorm0 = K0->GetMaxAbsEntry();

  double fNorm0 = 0.0;
  for(int j=0; j<r; j++)
    if (fabs(internalForce0[j]) > fNorm0)
      fNorm0 = fabs(internalForce0[j]);

  double eps = 1.0;
  do
  {
    for(int j=0; j<r; j++)
    {
      dqeps[j] = eps * dq[j];
      q1[j] = q[j] + dqeps[j];
    }

    GetForceAndMatrix(q1, internalForce1, K1);

    // internalForce1 - internalForce0 - K0 * dqeps should be O(eps^2)
    K0->MultiplyVector(dqeps, testV);

    for(int j=0; j<r; j++)
      testV[j] = internalForce1[j] - internalForce0[j] - testV[j];

    // compute largest abs value entry in testV
    double maxEntry = 0.0;
    for(int j=0; j<r; j++)
      if (fabs(testV[j]) > maxEntry)
        maxEntry = fabs(testV[j]);

    double dfNorm = 0.0;
    for(int j=0; j<r; j++)
      if (fabs(internalForce1[j] - internalForce0[j]) > dfNorm)
        dfNorm = fabs(internalForce1[j] - internalForce0[j]);
 
    printf("eps=%G: maxEntry=%G maxEntry/eps^2=%G dfNorm=%G f0Norm=%G K0Norm=%G\n", eps, maxEntry, maxEntry / eps / eps, dfNorm, fNorm0, KNorm0);
    eps *= 0.1;
  }
  while (eps > 1e-15); 

  delete(K0);

  free(testV);
  free(internalForce1);
  free(internalForce0);
  free(dqeps);
  free(q1);
}
