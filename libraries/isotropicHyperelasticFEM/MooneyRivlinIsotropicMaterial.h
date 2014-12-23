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

#ifndef _MOONEYRIVLINISOTROPICMATERIAL_H_
#define _MOONEYRIVLINISOTROPICMATERIAL_H_

#include "isotropicMaterialWithCompressionResistance.h"
#include "tetMesh.h"

/*
  The implemented compressible Mooney-Rivlin material has the following energy function:

  energy = 0.5 * mu01 * ((Ic^2 - IIc) / pow(IIIc, 2/3) - 6) +
                 mu10 * (Ic / pow(IIIc, 1/3) - 3) +
                 v1 * (sqrt(IIIc) - 1)^2

  This is a standard adaption of the material presented in:

  R. S. Rivlin and D. W. Saunders: 
  Large elastic deformations of isotropic materials VII. 
  Experiments on the deformation of rubber. 
  Philosophical Transactions Royal Society London Series A, 
  243(865), 1951, pp. 251-288.

  See the Wikipedia page, "Polynomial_(hyperelastic_model)".
  A snapshot of the Wikipedia page is included in this folder.
*/

class MooneyRivlinIsotropicMaterial : public IsotropicMaterialWithCompressionResistance
{
public:
  MooneyRivlinIsotropicMaterial(TetMesh * tetMesh, int enableCompressionResistance=0, double compressionResistance=0.0);
  virtual ~MooneyRivlinIsotropicMaterial();

  virtual double ComputeEnergy(int elementIndex, double * invariants);
  virtual void ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient); // invariants and gradient are 3-vectors
  virtual void ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian); // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).

protected:
  double * mu01_;
  double * mu10_;
  double * v1_;

  double compressionResistance;
  double * EdivNuFactor;
  virtual double GetCompressionResistanceFactor(int elementIndex);
};

#endif

