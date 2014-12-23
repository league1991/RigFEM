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
#include "neoHookeanIsotropicMaterial.h"
#include "volumetricMeshENuMaterial.h"

NeoHookeanIsotropicMaterial::NeoHookeanIsotropicMaterial(TetMesh * tetMesh, int enableCompressionResistance_, double compressionResistance_) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance_), compressionResistance(compressionResistance_)
{
  int numElements = tetMesh->getNumElements();
  lambdaLame = (double*) malloc (sizeof(double) * numElements);
  muLame = (double*) malloc (sizeof(double) * numElements);

  if (enableCompressionResistance)
    EdivNuFactor = (double*) malloc (sizeof(double) * numElements);
  else
    EdivNuFactor = NULL;

  for(int el=0; el<numElements; el++)
  {
    VolumetricMesh::Material * material = tetMesh->getElementMaterial(el);
    VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
    if (eNuMaterial == NULL)
    {
      printf("Error: NeoHookeanIsotropicMaterial: mesh does not consist of E, nu materials.\n");
      throw 1;
    }

    lambdaLame[el] = eNuMaterial->getLambda();
    muLame[el] = eNuMaterial->getMu();

    if (enableCompressionResistance)
    {
      EdivNuFactor[el] = compressionResistance * eNuMaterial->getE() / (1.0 - 2.0 * eNuMaterial->getNu());
      //printf("Setting EdivNuFactor[%d]=%G\n", el, EdivNuFactor[el]);
    }
  }
}

NeoHookeanIsotropicMaterial::~NeoHookeanIsotropicMaterial()
{
  free(EdivNuFactor);
  free(lambdaLame);
  free(muLame);
}

double NeoHookeanIsotropicMaterial::ComputeEnergy(int elementIndex, double * invariants)
{
  double IC = invariants[0];
  double IIIC = invariants[2];
  double J = sqrt(IIIC); 
  double logJ = log(J);
  // Note: computation of J and logJ will fail for an inverted element.
  // The IsotropicHyperelasticFEM class will prevent inversions (assuming proper
  // threshold was set), so normally this is not an issue.

  double energy = 0.5 * muLame[elementIndex] * (IC - 3.0) - muLame[elementIndex] * logJ + 0.5 * lambdaLame[elementIndex] * logJ * logJ;

  AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

  return energy;
}

void NeoHookeanIsotropicMaterial::ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient) // invariants and gradient are 3-vectors
{
  //printf("Entered NeoHookeanIsotropicMaterial::ComputeEnergyGradient\n");

  double IIIC = invariants[2];
  gradient[0] = 0.5 * muLame[elementIndex];
  gradient[1] = 0.0;
  gradient[2] = (-0.5 * muLame[elementIndex] + 0.25 * lambdaLame[elementIndex] * log(IIIC)) / IIIC;

  AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

void NeoHookeanIsotropicMaterial::ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
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
  hessian[5] = (0.25 * lambdaLame[elementIndex] + 0.5 * muLame[elementIndex] - 0.25 * lambdaLame[elementIndex] * log(IIIC)) / (IIIC * IIIC);

  AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}

double NeoHookeanIsotropicMaterial::GetCompressionResistanceFactor(int elementIndex)
{
  return EdivNuFactor[elementIndex];
}

