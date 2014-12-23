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
#include <stdio.h>
#include "reducedForceModel.h"

void ReducedForceModel::GetForceAndMatrix(double * u, double * internalForces, double * tangentStiffnessMatrix)
{
  GetInternalForce(u, internalForces);
  GetTangentStiffnessMatrix(u, tangentStiffnessMatrix);
}

void ReducedForceModel::TestStiffnessMatrix(int numTests, double qAmplitude)
{
  double * q = (double*) malloc (sizeof(double) * r);
  double * dq = (double*) malloc (sizeof(double) * r);
  double * q1 = (double*) malloc (sizeof(double) * r);
  double * f = (double*) malloc (sizeof(double) * r);
  double * df = (double*) malloc (sizeof(double) * r);
  double * f1 = (double*) malloc (sizeof(double) * r);
  double * K = (double*) malloc (sizeof(double) * r * r);
  for(int i=0; i<numTests; i++)
  {
    for(int j=0; j<r; j++)
    {
      q[j] = qAmplitude * (-1.0 + 2.0 * rand() / RAND_MAX);
      dq[j] = qAmplitude * (-1.0 + 2.0 * rand() / RAND_MAX);
    }

    GetInternalForce(q, f);
    GetTangentStiffnessMatrix(q, K);
    double eps = 1.0;
    do
    {
      for(int j=0; j<r; j++)
        q1[j] = q[j] + eps * dq[j];

      GetInternalForce(q1, f1);

      for(int j=0; j<r; j++)
      {
        df[j] = 0;
        for(int k=0; k<r; k++)
          df[j] += K[r*k+j] * dq[k];
        df[j] *= eps;
      }

      // f + eps * K * dq - f1 should be O(eps^2)
      double maxEntry = 0.0;
      for(int j=0; j<r; j++)
      {
        double entry = fabs( (f1[j] - f[j] - df[j]) / eps );
        if (entry > maxEntry)
          maxEntry = entry;
      }

      printf("Stiffness matrix trial %d: eps=%G: maxEntryntry=%G\n", i, eps, maxEntry);
      eps *= 0.5;
    }
    while (eps > 1e-15);
  }

  free(K);
  free(dq);
  free(q1);
  free(q);
  free(f);
  free(df);
  free(f1);
}

