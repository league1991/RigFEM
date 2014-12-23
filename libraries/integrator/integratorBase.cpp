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

#include "integratorBase.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

IntegratorBase::IntegratorBase(int r, double timestep, double dampingMassCoef, double dampingStiffnessCoef)
{
  this->r = r;
  this->timestep = timestep;

  this->dampingMassCoef = dampingMassCoef;
  this->dampingStiffnessCoef = dampingStiffnessCoef;

  internalForceScalingFactor = 1.0;

  q = (double*) malloc(sizeof(double) * r);
  qvel = (double*) malloc(sizeof(double) * r);
  qaccel = (double*) malloc(sizeof(double) * r);
  q_1 = (double*) malloc(sizeof(double) * r);
  qvel_1 = (double*) malloc(sizeof(double) * r);
  qaccel_1 = (double*) malloc(sizeof(double) * r);
  externalForces = (double*) malloc(sizeof(double) * r);
  internalForces = (double*) malloc(sizeof(double) * r);

  qresidual = (double*) malloc(sizeof(double) * r);
  qdelta = (double*) malloc(sizeof(double) * r);
  buffer = (double*) malloc(sizeof(double) * r);

  ResetToRest();

  memset(externalForces,0,sizeof(double) * r);
  memset(internalForces,0,sizeof(double) * r);
}

IntegratorBase::~IntegratorBase()
{
  free(q);
  free(qvel);
  free(qaccel);
  free(q_1);
  free(qvel_1);
  free(qaccel_1);

  free(externalForces);
  free(internalForces);

  free(qresidual);
  free(qdelta);
  free(buffer);
}

void IntegratorBase::SetExternalForces(double * externalForces)
{
  memcpy(this->externalForces,externalForces,sizeof(double)*r);
}

void IntegratorBase::AddExternalForces(double * externalForces)
{
  for(int i=0; i<r; i++)
    this->externalForces[i] += externalForces[i];
}

void IntegratorBase::SetExternalForcesToZero()
{
  memset(externalForces,0,sizeof(double)*r);
}

void IntegratorBase::GetExternalForces(double * externalForces_copy)
{
  memcpy(externalForces_copy,externalForces,sizeof(double)*r);
}

void IntegratorBase::SetqState(const double * q, const double * qvel, const double * qaccel)
{
  memcpy(this->q,q,sizeof(double)*r);
  if (qvel != NULL)
    memcpy(this->qvel,qvel,sizeof(double)*r);
  if (qaccel != NULL)
    memcpy(this->qaccel,qaccel,sizeof(double)*r);
}

void IntegratorBase::GetqState(double * q, double * qvel, double * qaccel)
{
  if (q != NULL)
    memcpy(q,this->q,sizeof(double)*r);
  if (qvel != NULL)
    memcpy(qvel,this->qvel,sizeof(double)*r);
  if (qaccel != NULL)
    memcpy(qaccel,this->qaccel,sizeof(double)*r);
}

void IntegratorBase::ResetToRest()
{
  memset(q,0,sizeof(double)*r);
  memset(qvel,0,sizeof(double)*r);
  memset(qaccel,0,sizeof(double)*r);
  memset(q_1,0,sizeof(double)*r);
  memset(qvel_1,0,sizeof(double)*r);
  memset(qaccel_1,0,sizeof(double)*r);
}

void IntegratorBase::ConstrainToSphere(double R2)
{
  double norm2 = 0;
  for(int i=0; i<r; i++)
    norm2 += q[i] * q[i];

  if (norm2 > R2)
  {
    double beta = sqrt(R2 / norm2);
    for(int dim=0; dim<r; dim++)
    {
      q[dim] *= beta;
      q_1[dim] = q[dim];
      qvel[dim] = 0;
      qvel_1[dim] = 0;
      qaccel[dim] = 0;
      qaccel_1[dim] = 0;
    }
  }
}

