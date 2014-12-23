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
#include "StVKIsotropicMaterial.h"
#include "volumetricMeshENuMaterial.h"

StVKIsotropicMaterial::StVKIsotropicMaterial(TetMesh * tetMesh, int enableCompressionResistance_, double compressionResistance_) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance_), compressionResistance(compressionResistance_)
{
  //printf("Entering StVKIsotropicMaterial::StVKIsotropicMaterial\n");
  //printf("Enable compression resistance is: %d. Compression resistance is: %G\n", enableCompressionResistance, compressionResistance);
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
      printf("Error: StVKIsotropicMaterial: mesh does not consist of E, nu materials.\n");
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

StVKIsotropicMaterial::~StVKIsotropicMaterial() 
{
  free(EdivNuFactor);
  free(muLame);
  free(lambdaLame);
}

double StVKIsotropicMaterial::ComputeEnergy(int elementIndex, double * invariants)
{
  double IC = invariants[0];
  double IIC = invariants[1];
  //double IIIC = invariants[2]; // not needed for StVK

  double energy = 0.125 * lambdaLame[elementIndex] * (IC - 3.0) * (IC - 3.0) + 0.25 * muLame[elementIndex] * (IIC - 2.0 * IC + 3.0);

  AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

/*
  if (enableInversionPrevention)
  {
    double IIIC = invariants[2]; 
    double J = sqrt(IIIC);

    if (J < 1)
    {
      double fac = (J - 1) * (J - 1) * (J - 1) / (6 * 6 * 6);
      energy += -EdivNuFactor[elementIndex] * fac / 12;
    }
  }  
*/

  return energy;
}

void StVKIsotropicMaterial::ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient) // invariants and gradient are 3-vectors
{
  //printf("Entered StVKIsotropicMaterial::ComputeEnergyGradient\n");

  double IC = invariants[0];
  gradient[0] = 0.25 * lambdaLame[elementIndex] * (IC - 3.0) - 0.5 * muLame[elementIndex];
  gradient[1] = 0.25 * muLame[elementIndex];
  gradient[2] = 0.0;

  AddCompressionResistanceGradient(elementIndex, invariants, gradient);

/*
  if (enableInversionPrevention)
  {
    double IIIC = invariants[2]; 
    double J = sqrt(IIIC);

    if (J < 1)
    {
      double fac = (J - 1) * (J - 1) / (6 * 6);
      gradient[2] += -EdivNuFactor[elementIndex] * fac / (8 * J);
    }
  }  
*/
}

void StVKIsotropicMaterial::ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
{
  // 11
  hessian[0] = 0.25 * lambdaLame[elementIndex];
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

/*
  if (enableInversionPrevention)
  {
    double IIIC = invariants[2]; 
    double J = sqrt(IIIC);

    if (J < 1)
      hessian[5] += EdivNuFactor[elementIndex] * (1 - J) * (1 + 11 * J) / (576 * J * J * J);
  }  
*/
}

double StVKIsotropicMaterial::GetCompressionResistanceFactor(int elementIndex)
{
  return EdivNuFactor[elementIndex];
}

