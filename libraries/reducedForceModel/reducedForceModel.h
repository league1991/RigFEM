/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "reducedForceModel" library , Copyright (C) 2007 CMU, 2009 MIT,       *
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
  Abstract class for a reduced deformable material model.
  Serves as a connecting class between reduced integrators and classes to calculate reduced internal forces and tangent stiffness matrices.
*/

#ifndef _REDUCEDFORCEMODELMODEL_H_
#define _REDUCEDFORCEMODELMODEL_H_

#include <stdlib.h>
#include <stdio.h>

class ReducedForceModel
{
public:
  ReducedForceModel() {}
  virtual ~ReducedForceModel() {}

  inline int Getr() { return r;}

  virtual void GetInternalForce(double * u, double * internalForces) = 0;
  virtual void GetTangentStiffnessMatrix(double * u, double * tangentStiffnessMatrix) = 0; 
  // sometimes computation time can be saved if we know that we will need both internal forces and tangent stiffness matrices:
  virtual void GetForceAndMatrix (double * u, double * internalForces, double * tangentStiffnessMatrix); 
  virtual void ResetToZero() {}
  virtual void Reset(double * q) {}

  void TestStiffnessMatrix(int numTrials, double qMagnitude = 1.0);

protected:
  int r;
};


#endif

