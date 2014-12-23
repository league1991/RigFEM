/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "elasticForceModel" library , Copyright (C) 2007 CMU, 2009 MIT,       *
 *                                                       2014 USC        *
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

/*
  Nonlinear reduced StVK force model
  
  The reduced force and reduced stiffness matrix are computed in the 
  high-dimensional space, and then projected:
  fq = U^T f(U*q)
*/

#ifndef _REDUCEDSUBSPACESTVKFORCEMODEL_H_
#define _REDUCEDSUBSPACESTVKFORCEMODEL_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "reducedForceModel.h"
#include "StVKStiffnessMatrix.h"
#include "StVKInternalForces.h"
#include "modalMatrix.h"
#include "sparseMatrix.h"

class ReducedSubspaceStVKForceModel : public ReducedForceModel
{
public:
  
  ReducedSubspaceStVKForceModel(StVKInternalForces * stVKInternalForces, 
                                StVKStiffnessMatrix * stVKStiffnessMatrix,
                                ModalMatrix * modalMatrix);
  virtual ~ReducedSubspaceStVKForceModel();

  virtual void GetInternalForce(double * q, double * internalForces); 
  virtual void GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix);
  virtual void GetForceAndMatrix(double * q, 
    double * internalForces, double * tangentStiffnessMatrix); 

protected:

  void GetTangentStiffnessMatrixHelper(double * tangentStiffnessMatrix);

  StVKInternalForces * stVKInternalForces;
  StVKStiffnessMatrix * stVKStiffnessMatrix;
  SparseMatrix * sparseMatrix;
  ModalMatrix * modalMatrix;
  int r;
  int n;
  double * U;
  double * u;
  double * bufferVector;
  double * bufferMatrix; 
};

#endif

