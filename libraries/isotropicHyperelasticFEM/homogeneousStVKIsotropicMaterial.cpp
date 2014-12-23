/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "isotropic hyperelastic FEM" library , Copyright (C) 2014 USC         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin                            *
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
#include "homogeneousStVKIsotropicMaterial.h"

HomogeneousStVKIsotropicMaterial::HomogeneousStVKIsotropicMaterial(double E_, double nu_, int enableCompressionResistance_, double compressionResistance_) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance_), compressionResistance(compressionResistance_)
{
  SetYoungModulusAndPoissonRatio(E_, nu_);
  EdivNuFactor = compressionResistance * E / (1.0 - 2.0 * nu);
}

HomogeneousStVKIsotropicMaterial::~HomogeneousStVKIsotropicMaterial() {}

void HomogeneousStVKIsotropicMaterial::SetYoungModulusAndPoissonRatio(double E_, double nu_)
{
  E = E_;
  nu = nu_;
  lambdaLame = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
  muLame = E / (2.0 * (1.0 + nu));
}

void HomogeneousStVKIsotropicMaterial::SetLameCoefficients(double lambda_, double mu_)
{
  lambdaLame = lambda_;
  muLame = mu_;
}

double HomogeneousStVKIsotropicMaterial::ComputeEnergy(int elementIndex, double * invariants)
{
  double IC = invariants[0];
  double IIC = invariants[1];
  //double IIIC = invariants[2]; // not needed for StVK

  double energy = 0.125 * lambdaLame * (IC - 3.0) * (IC - 3.0) + 0.25 * muLame * (IIC - 2.0 * IC + 3.0);

  AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

  return energy;
}

void HomogeneousStVKIsotropicMaterial::ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient) // invariants and gradient are 3-vectors
{
  double IC = invariants[0];
  gradient[0] = 0.25 * lambdaLame * (IC - 3.0) - 0.5 * muLame;
  gradient[1] = 0.25 * muLame;
  gradient[2] = 0.0;

  AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

void HomogeneousStVKIsotropicMaterial::ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
{
  // 11
  hessian[0] = 0.25 * lambdaLame;
  // 12
  hessian[1] = 0.0;
  // 13
  hessian[2] = 0.0;
  // 22
  hessian[3] = 0.0;
  // 23
  hessian[4] = 0.0;
  // 33
  hessian[5] = 0.0;

  AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}

double HomogeneousStVKIsotropicMaterial::GetCompressionResistanceFactor(int elementIndex)
{
  return EdivNuFactor;
}

