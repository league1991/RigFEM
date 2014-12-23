/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "reducedStVK" library , Copyright (C) 2007 CMU, 2009 MIT              *
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

#ifndef _STVKREDUCEDHESSIANTENSOR_H_
#define _STVKREDUCEDHESSIANTENSOR_H_

#include "StVKReducedStiffnessMatrix.h"

/*
  This class computes the reduced Hessian tensor (gradient of the reduced stiffness matrix),
  for a reduced St. Venant Kirchhoff deformable model.
  The tensor is rank 3; its dimensions are r x r x r. 
  Each element of the tensor is a linear polynomial in the reduced coordinates q.

  See also StVKReducedInternalForces.h .
*/

class StVKReducedHessianTensor
{
public:

  // creates the linear coefficients of the tensor, using those from the reduced stiffness matrix class
  StVKReducedHessianTensor(StVKReducedStiffnessMatrix * stVKReducedStiffnessMatrix);

  // load the coefficients from file
  StVKReducedHessianTensor(const char * filename);

  ~StVKReducedHessianTensor();

  // saves coefficients out, in binary format
  int Save(const char * filename);

  // make a buffer that you can then pass it to the "Evaluate" routine
  void MakeRoomForTensor(double ** Hq);

  // evaluates the Hessian tensor matrix for the given configuration q, result is written into Hq
  void Evaluate(double * q, double * Hq);

  // computes A = Hq : q 
  // (A is a r x r matrix; it will be symmetric (however all r*r entries are returned, not just upper triangle))
  static void ContractWithVector(int r, double * Hq, double * q, double * A);

  inline int Getr() { return r; }
  // print the Hessian tensor to standard output, in Mathematica format
  void PrintTensor();

  // === advanced routines below here ===

  // evaluates only entries [start,...end), in the upper-triangular order 
  //   (0 <= start <= end <= quadraticSize)
  void EvaluateSubset(double * q, int start, int end, double * Hq);

  // approximate quantities using a local hessian model
  // output fields: dfq, dRq
  void ApproximateReducedForceAndStiffnessMatrix
    (double * baseHessianTensor, double * baseStiffnessMatrix, 
     double * dq, double * dfq, double * dK);

  inline double freeCoef(int output, int input, int deriv) { return freeCoef_[freeCoefPos(output, input, deriv)]; }
  inline double linearCoef(int output, int input, int deriv, int i) { return linearCoef_[linearCoefPos(output, input, deriv, i)]; }

  // makes shallow copies of all pointers, except those initialized by InitBuffers
  // use this if you want to Evaluate two or more identical models (i.e., two copies of an object) in parallel (to ensure thread safety)
  // you do not need to use this if you are Evaluating a single model in parallel (e.g., using the MT derived class)
  StVKReducedHessianTensor * ShallowClone();

protected:

/*
  The Hessian is internally stored as follows:
  rrrrrrr
   rrrrrr
    rrrrr
     rrrr
      rrr
       rr
        r

  Each "r" is a vector of length r.

  Hessian coefficients are internally stored in the following format:
  
  free terms:
  rrrrrrr
   rrrrrr
    rrrrr
     rrrr
      rrr
       rr
        r

  linear terms:
  qqqqqqq
   qqqqqq
    qqqqq
     qqqq
      qqq
       qq
        q

  Each "q" is a r x r matrix of vectors of length r:
  rrrrrrr
  rrrrrrr
  rrrrrrr
  rrrrrrr
  rrrrrrr
  rrrrrrr
  rrrrrrr
*/
  
  int r,r2;

  double * freeCoef_;
  double * linearCoef_;

  double * buffer1, * buffer2;

  int linearSize, quadraticSize;

  inline int tensorPosition(int output, int input, int deriv);

  // assumes input >= output, deriv is not constrained
  inline int freeCoefPos(int output, int input, int deriv) { return tensorPosition(output,input,deriv); }
  // assumes input >= output, deriv and i are not constrained
  inline int linearCoefPos(int output, int input, int deriv, int i) {return tensorPosition(output,input,deriv) * linearSize + i;}

  void InitBuffers();
  void FreeBuffers();

  int shallowCopy;
};

inline int StVKReducedHessianTensor::tensorPosition(int i, int j, int k)
{ 
  return (i * r - i * (i+1) / 2 + j) * r + k;
}

#endif

