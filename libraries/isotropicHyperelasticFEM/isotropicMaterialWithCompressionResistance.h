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

#ifndef _ISOTROPICMATERIALWITHCOMPRESSIONRESISTANCE_H_
#define _ISOTROPICMATERIALWITHCOMPRESSIONRESISTANCE_H_

#include "isotropicMaterial.h"

class IsotropicMaterialWithCompressionResistance : public IsotropicMaterial
{
public:
  IsotropicMaterialWithCompressionResistance(int enableCompressionResistance=0);
  virtual ~IsotropicMaterialWithCompressionResistance();

  virtual double ComputeEnergy(int elementIndex, double * invariants)=0;
  virtual void ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient)=0; // invariants and gradient are 3-vectors
  virtual void ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian)=0; // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).

protected:
  int enableCompressionResistance;

  virtual double GetCompressionResistanceFactor(int elementIndex);
  void AddCompressionResistanceEnergy(int elementIndex, double * invariants, double * energy);
  void AddCompressionResistanceGradient(int elementIndex, double * invariants, double * gradient);
  void AddCompressionResistanceHessian(int elementIndex, double * invariants, double * hessian); 
};

#endif

