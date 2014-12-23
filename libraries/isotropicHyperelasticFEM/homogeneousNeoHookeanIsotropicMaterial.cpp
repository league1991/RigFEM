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
#include "homogeneousNeoHookeanIsotropicMaterial.h"

HomogeneousNeoHookeanIsotropicMaterial::HomogeneousNeoHookeanIsotropicMaterial(double E_, double nu_, int enableCompressionResistance_, double compressionResistance_) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance_), compressionResistance(compressionResistance_)
{
  SetYoungModulusAndPoissonRatio(E_, nu_);
  EdivNuFactor = compressionResistance * E / (1.0 - 2.0 * nu);
}

HomogeneousNeoHookeanIsotropicMaterial::~HomogeneousNeoHookeanIsotropicMaterial() {}

void HomogeneousNeoHookeanIsotropicMaterial::SetYoungModulusAndPoissonRatio(double E_, double nu_)
{
  E = E_;
  nu = nu_;
  lambdaLame = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
  muLame = E / (2.0 * (1.0 + nu));
}

void HomogeneousNeoHookeanIsotropicMaterial::SetLameCoefficients(double lambda_, double mu_)
{
  lambdaLame = lambda_;
  muLame = mu_;
}

double HomogeneousNeoHookeanIsotropicMaterial::ComputeEnergy(int elementIndex, double * invariants)
{
  double IC = invariants[0];
  double IIIC = invariants[2];
  double J = sqrt(IIIC); 
  double logJ = log(J);
  // Note: computation of J and logJ will fail for an inverted element.
  // The IsotropicHyperelasticFEM class will prevent inversions (assuming proper
  // threshold was set), so normally this is not an issue.

  double energy = 0.5 * muLame * (IC - 3.0) - muLame * logJ + 0.5 * lambdaLame * logJ * logJ;

  AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

  return energy;
}

void HomogeneousNeoHookeanIsotropicMaterial::ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient) // invariants and gradient are 3-vectors
{
  double IIIC = invariants[2];

  gradient[0] = 0.5 * muLame;
  gradient[1] = 0.0;
  gradient[2] = (-0.5 * muLame + 0.25 * lambdaLame * log(IIIC)) / IIIC;

  AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

void HomogeneousNeoHookeanIsotropicMaterial::ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
{
  double IIIC = invariants[2];
  // 11
  hessian[0] = 0.0;
  // 12
  hessian[1] = 0.0;
  // 13
  hessian[2] = 0.0;
  // 22
  hessian[3] = 0.0;
  // 23
  hessian[4] = 0.0;
  // 33
  hessian[5] = (0.25 * lambdaLame + 0.5 * muLame - 0.25 * lambdaLame * log(IIIC)) / (IIIC * IIIC);

  AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}

double HomogeneousNeoHookeanIsotropicMaterial::GetCompressionResistanceFactor(int elementIndex)
{
  return EdivNuFactor;
}

