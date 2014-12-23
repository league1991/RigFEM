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
#include <math.h>
#include "lapack-headers.h"
#include "matrix.h"
#include "matrixLAPACK.h"
#include "performanceCounter.h"
#include "implicitNewmarkDense.h"
#include "implicitNewmarkDenseMulti1D.h"
#include "IPIVC.h"

ImplicitNewmarkDenseMulti1D::ImplicitNewmarkDenseMulti1D(int r_, double timestep_, double * massMatrix_, double * tangentStiffnessMatrix_, double dampingMassCoef_, double dampingStiffnessCoef_, double NewmarkBeta_, double NewmarkGamma_):
  IntegratorBaseDense(r_, timestep_, dampingMassCoef_, dampingStiffnessCoef_),
  IntegratorMulti1D(r_, timestep_, massMatrix_, tangentStiffnessMatrix_, dampingMassCoef_, dampingStiffnessCoef_),
  ImplicitNewmarkDense(r_, timestep_, dampingMassCoef_, dampingStiffnessCoef_, 0.25, 0.5)
{
  coef_deltaZ = (double *) malloc (sizeof(double) * r);
  coef_zvel = (double *) malloc (sizeof(double) * r);
  coef_zaccel = (double *) malloc (sizeof(double) * r);

  UpdateCoefs();
}

ImplicitNewmarkDenseMulti1D::~ImplicitNewmarkDenseMulti1D()
{
  free(coef_deltaZ);
  free(coef_zvel);
  free(coef_zaccel);
}

int ImplicitNewmarkDenseMulti1D::DoTimestep()
{
  // rotate the external forces: rotatedForces = Q^T * externalForces
  bool useRotationTranspose = true;
  RotateVector(Q, externalForces, rotatedForces, useRotationTranspose);

  for(int dim=0; dim<r; dim++)
  {
    // evaluate rhs = (alpha2 - damping * alpha5) * z_vel + (alpha3 - damping * alpha6) * z_accel - stiffness * z + f_ext
    double rhs = coef_zvel[dim] * zvel[dim] + coef_zaccel[dim] * zaccel[dim] - tangentStiffnessMatrix[dim] * z[dim] + rotatedForces[dim];

    // solve for qDelta
    qdelta[dim] = rhs / coef_deltaZ[dim];
  }

  // update z, zvel, and zaccel
  for(int dim=0; dim<r; dim++)
  {
    z[dim] += qdelta[dim];
    double zvelTemp = zvel[dim];
    double zaccelTemp = zaccel[dim];
    zvel[dim] = alpha4 * qdelta[dim] + alpha5 * zvelTemp + alpha6 * zaccelTemp;
    zaccel[dim] = alpha1 * qdelta[dim] - alpha2 * zvelTemp - alpha3 * zaccelTemp;
  }

  // update states by rotating z, zvel, zaccel with Q
  useRotationTranspose = false;
  RotateVector(Q, z, q, useRotationTranspose); // q = Q z
  RotateVector(Q, zvel, qvel, useRotationTranspose); // qvel = Q zvel
  RotateVector(Q, zaccel, qaccel, useRotationTranspose); // qaccel = Q zaccel

  return 0;
}

void ImplicitNewmarkDenseMulti1D::SetTimestep(double timestep_)
{
  ImplicitNewmarkDense::SetTimestep(timestep_);
  UpdateCoefs();
}

// this function should be always called after alphas are updated
void ImplicitNewmarkDenseMulti1D::UpdateCoefs()
{
  for(int dim=0; dim<r; dim++)
  {
    dampingMatrix[dim] = dampingMassCoef + dampingStiffnessCoef * tangentStiffnessMatrix[dim];
    coef_deltaZ[dim] = alpha1 + dampingMatrix[dim] * alpha4 + tangentStiffnessMatrix[dim];
    coef_zvel[dim] = alpha2 - dampingMatrix[dim] * alpha5;
    coef_zaccel[dim] = alpha3 - dampingMatrix[dim] * alpha6;
  }
}

void ImplicitNewmarkDenseMulti1D::SetNewmarkBeta(double NewmarkBeta_)
{
  NewmarkBeta = NewmarkBeta_; 
  UpdateAlphas(); 
  UpdateCoefs();
}

void ImplicitNewmarkDenseMulti1D::SetNewmarkGamma(double NewmarkGamma_)
{
  NewmarkGamma = NewmarkGamma_; 
  UpdateAlphas();
  UpdateCoefs();
}

void ImplicitNewmarkDenseMulti1D::SetInternalForceScalingFactor(double internalForceScalingFactor_)
{
  internalForceScalingFactor = internalForceScalingFactor_;
  for(int dim=0; dim<r; dim++)
    tangentStiffnessMatrix[dim] = internalForceScalingFactor * tangentStiffnessMatrixOriginal[dim];
  UpdateCoefs();
}

void ImplicitNewmarkDenseMulti1D::SetDampingMassCoef(double dampingMassCoef_)
{ 
  dampingMassCoef = dampingMassCoef_;
  UpdateCoefs();
}

void ImplicitNewmarkDenseMulti1D::SetDampingStiffnessCoef(double dampingStiffnessCoef_)
{
  dampingStiffnessCoef = dampingStiffnessCoef_;
  UpdateCoefs();
}

