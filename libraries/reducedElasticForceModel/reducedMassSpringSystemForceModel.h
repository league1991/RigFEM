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
  Force model for a reduced mass spring system.
  The reduced force and reduced stiffness matrix are computed in the 
  high-dimensional space, and then projected:
  fq = U^T f(U*q)
*/

#ifndef _REDUCEDMASSSPRINGSYSTEMFORCEMODEL_H_
#define _REDUCEDMASSSPRINGSYSTEMFORCEMODEL_H_

#include <stdlib.h>
#include "modalMatrix.h"
#include "sparseMatrix.h"
#include "massSpringSystem.h"
#include "reducedForceModel.h"

class ReducedMassSpringSystemForceModel : public virtual ReducedForceModel
{
public:
  ReducedMassSpringSystemForceModel(MassSpringSystem * massSpringSystem, ModalMatrix * modalMatrix);

  virtual ~ReducedMassSpringSystemForceModel();
  virtual void GetInternalForce(double * q, double * internalForces); 
  virtual void GetTangentStiffnessMatrix(double * q, double * tangentStiffnessMatrix);

  virtual void GetForceAndMatrix(double * q,
    double * internalForces, double * tangentStiffnessMatrix);

protected:
  void GetTangentStiffnessMatrixHelper(double * tangentStiffnessMatrix);

  MassSpringSystem * massSpringSystem;
  ModalMatrix * modalMatrix;
  SparseMatrix * sparseMatrix;

  int n;
  double * U;
  double * u;
  double * bufferVector;
  double * bufferMatrix;

  // this macro is normally defined in the makefile (header)
  #if USE_MKL_SPARSE_BLAS
    double * csr_values;
    int * csr_columns;
    int * csr_pointerB;
    int * csr_pointerE;
  #endif
};

#endif

