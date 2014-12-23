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
#include "homogeneousMooneyRivlinIsotropicMaterial.h"

HomogeneousMooneyRivlinIsotropicMaterial::HomogeneousMooneyRivlinIsotropicMaterial(double mu01_, double mu10_, double v1_, int enableCompressionResistance_, double compressionResistance_) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance_), mu01(mu01_), mu10(mu10_), v1(v1_), compressionResistance(compressionResistance_) 
{
  // invert the following formulas to compute "pseudo" E and nu that correspond to this material
  // (only needed for compression resistance)

  //K = Ea / ( 3 * ( 1 - 2 * nua ) );
  //G = Ea / ( 2 * ( 1 + nua ) );
  //mu10 = G / ( 2 * ( 1 + mur ) );
  //mu01 = mur * mu10;
  //v1 = K / 2;

  double K = 2.0 * v1;
  double muRatio = mu01 / mu10;
  double G = mu10 * ( 2 * ( 1 + muRatio ) );

  double E = 9 * K * G / (3 * K + G);
  double nu = (3 * K - 2 * G) / (2 * (3 * K + G));

  EdivNuFactor = compressionResistance * E / (1.0 - 2.0 * nu);
}

HomogeneousMooneyRivlinIsotropicMaterial::~HomogeneousMooneyRivlinIsotropicMaterial() {}

double HomogeneousMooneyRivlinIsotropicMaterial::ComputeEnergy(int elementIndex, double * invariants)
{
  double Ic = invariants[0];
  double IIc = invariants[1];
  double IIIc = invariants[2];
  double energy = 0.5 * (-6.0 + (Ic * Ic - IIc) / pow(IIIc, 2.0 / 3.0)) * mu01 + 
                  (-3.0 + Ic / pow(IIIc, 1.0 / 3.0)) * mu10 + 
                  pow(-1.0 + sqrt(IIIc), 2.0) * v1;

  AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

  return energy;
}

void HomogeneousMooneyRivlinIsotropicMaterial::ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient) // invariants and gradient are 3-vectors
{
  double Ic = invariants[0];
  double IIc = invariants[1];
  double IIIc = invariants[2];
  gradient[0] = (Ic * mu01) / pow(IIIc, 2.0 / 3.0) + 
    mu10 / pow(IIIc, 1.0 / 3.0);
  gradient[1] = (-0.5 * mu01) / pow(IIIc, 2.0 / 3.0);
  gradient[2] = (-1.0 / 3.0 * (Ic * Ic - IIc) * mu01) / pow(IIIc, 5.0 / 3.0) - 
    (1.0 / 3.0 * Ic * mu10) / pow(IIIc, 4.0 / 3.0) + 
    ((-1.0 + sqrt(IIIc)) * v1) / sqrt(IIIc);

  AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

void HomogeneousMooneyRivlinIsotropicMaterial::ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
{
  double Ic = invariants[0];
  double IIc = invariants[1];
  double IIIc = invariants[2];

  // 11
  hessian[0] = mu01 / pow(IIIc, 2.0 / 3.0);
  // 12
  hessian[1] = 0.0;
  // 13
  hessian[2] = (-2.0 / 3.0) * Ic * mu01 / pow(IIIc, 5.0 / 3.0) - 
               mu10 / (3.0 * pow(IIIc, 4.0 / 3.0));
  // 22
  hessian[3] = 0.0;
  // 23
  hessian[4] = mu01 / (3.0 * pow(IIIc, 5.0 / 3.0));
  // 33
  hessian[5] = (5.0 / 9.0) * (Ic * Ic - IIc) * mu01 / pow(IIIc, 8.0 / 3.0) + 
               ((4.0 / 9.0) * Ic * mu10) / pow(IIIc, 7.0 / 3.0) - 
               (-1.0 + sqrt(IIIc)) * v1 / (2.0 * pow(IIIc, 1.5)) + 
               v1 / (2.0 * IIIc);

  AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}

double HomogeneousMooneyRivlinIsotropicMaterial::GetCompressionResistanceFactor(int elementIndex)
{
  return EdivNuFactor;
}

