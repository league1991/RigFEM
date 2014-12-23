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
#include "isotropicMaterialWithCompressionResistance.h"

IsotropicMaterialWithCompressionResistance::IsotropicMaterialWithCompressionResistance(int enableCompressionResistance_) : IsotropicMaterial(), enableCompressionResistance(enableCompressionResistance_)
{
}

IsotropicMaterialWithCompressionResistance::~IsotropicMaterialWithCompressionResistance() {}

void IsotropicMaterialWithCompressionResistance::AddCompressionResistanceEnergy(int elementIndex, double * invariants, double * energy)
{
  if (enableCompressionResistance)
  {
    double IIIC = invariants[2];
    double J = sqrt(IIIC);

    if (J < 1)
    {
      double compressionResistanceFactor = GetCompressionResistanceFactor(elementIndex);
      *energy += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) * (J - 1.0) / 2592.0;
    }
  }
}

void IsotropicMaterialWithCompressionResistance::AddCompressionResistanceGradient(int elementIndex, double * invariants, double * gradient)
{
  if (enableCompressionResistance)
  {
    double IIIC = invariants[2];
    double J = sqrt(IIIC);

    if (J < 1)
    {
      double compressionResistanceFactor = GetCompressionResistanceFactor(elementIndex);
      gradient[2] += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
    }
  }
}

void IsotropicMaterialWithCompressionResistance::AddCompressionResistanceHessian(int elementIndex, double * invariants, double * hessian)
{
  if (enableCompressionResistance)
  {
    double IIIC = invariants[2]; 
    double J = sqrt(IIIC);

    if (J < 1.0)
    {
      double compressionResistanceFactor = GetCompressionResistanceFactor(elementIndex);
      hessian[5] += compressionResistanceFactor * (1.0 - J) * (1.0 + J) / (3456.0 * J * J * J);
    }
  }
}

double IsotropicMaterialWithCompressionResistance::GetCompressionResistanceFactor(int elementIndex)
{
  return 1.0; // generic
}

