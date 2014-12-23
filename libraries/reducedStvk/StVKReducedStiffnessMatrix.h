/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "reducedStvk" library , Copyright (C) 2007 CMU, 2009 MIT              *
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
  This class computes the reduced tangent stiffness matrices for a reduced 
  St. Venant Kirchhoff deformable model. Each element of the stiffness matrix 
  is a quadratic polynomial in the reduced coordinates q.
  The stiffness matrix is a symmetric r x r matrix. Reduced coordinates 
  are a vector of length r, where r is the dimension of the reduced basis.

  The class serves two functions:
  1. Precompute the coefficients of the quadratic polynomials, given an 
     instance of the reduced internal forces class
  2. Given the coefficients, evaluate the reduced tangent stiffness matrix quickly. 
     The evaluation routine uses BLAS Level 2 and 3 routines. 
     The evaluation routine has been optimized extensively and is very efficient.

  See also StVKReducedInternalForces.h .
*/


#ifndef _STVKREDUCEDSTIFFNESSMATRIX_H_
#define _STVKREDUCEDSTIFFNESSMATRIX_H_

#include "StVKReducedInternalForces.h"

class StVKReducedStiffnessMatrix
{
public:
  // creates the reduced stiffness matrix coefficients, using a reduced internal forces class
  // note: this computation is very fast; it is much faster than computing the reduced internal force polynomials
  // therefore, loading and saving these coefficients is rarely necessary in practice
  StVKReducedStiffnessMatrix(StVKReducedInternalForces * stVKReducedInternalForces, int verbose=1);

  // load the coefficients from a file
  StVKReducedStiffnessMatrix(const char * filename);

  ~StVKReducedStiffnessMatrix();

  // saves coefficients (in binary format; the convention is to use the .sti file extension)
  int Save(const char * filename);

  // evaluates the stiffness matrix for the given configuration q, result is written into Kq (must be a pre-allocated r x r matrix)
  // Kq will be symmetric; the routine returns all the r*r entries of the matrix, as opposed to just the upper triangle
  void Evaluate(double * q, double * Kq);
  // evaluates the matrix assuming a purely linear model; the result is the reduced stiffness matrix in the rest configuration
  // note: parameter q does not affect the output in this case
  void EvaluateLinear(double * q, double * Kq); 

  inline int Getr() { return r; }

  // prints the matrix out to standard output, in Mathematica format
  void PrintMatrix();

  // ==== all routines below here are low-level research routines (expected to be rarely used) ===
  
  // query for coefficients:
  // assumes input>=output
  inline double freeCoef(int output, int input) { return freeCoef_[freeCoefPos(output,input)]; }
  // assumes input>=output
  inline double linearCoef(int output, int input, int i) { return linearCoef_[linearCoefPos(output,input,i)]; }
  // assumes (input>=output) && (j>=i)
  inline double quadraticCoef(int output, int input, int i, int j) { return quadraticCoef_[quadraticCoefPos(output,input,i,j)]; }

  // evaluates only entries [start,...end), in the upper-triangular order 
  //   (0 <= start <= end <= quadraticSize)
  void EvaluateSubset(double * q, int start, int end, double * Kq);

  void UseSingleThread(int useSingleThread);

  // makes shallow copies of all pointers, except those initialized by InitBuffers
  // use this if you want to Evaluate two or more identical models (i.e., two copies of an object) in parallel (to ensure thread safety)
  // you do not need to use this if you are Evaluating a single model in parallel (e.g., using the MT derived class)
  StVKReducedStiffnessMatrix * ShallowClone();

protected:

  // The stiffness matrix coefficients are stored as follows:

  // free terms:
  // 1111111
  //  111111
  //   11111
  //    1111
  //     111
  //      11
  //       1

  // linear terms:
  // rrrrrrr
  //  rrrrrr
  //   rrrrr
  //    rrrr
  //     rrr
  //      rr
  //       r

  // Each "r" is a vector of length r.

  // quadratic terms:
  // qqqqqqq
  //  qqqqqq
  //   qqqqq
  //    qqqq
  //     qqq
  //      qq
  //       q

  // Each "q" is an upper triangle of vectors of length r:
  // rrrrrrr
  //  rrrrrr
  //   rrrrr
  //    rrrr
  //     rrr
  //      rr
  //       r

  int r,r2;

  double * freeCoef_;
  double * linearCoef_;
  double * quadraticCoef_;

  int linearSize;
  int quadraticSize;

  inline int matrixPosition(int output, int input) { return output * r - output * (output+1) / 2 + input; }

  inline int freeCoefPos(int output, int input) { return matrixPosition(output,input); }
  inline int linearCoefPos(int output, int input, int i) {return matrixPosition(output,input) * linearSize + i;}
  inline int quadraticCoefPos(int output, int input, int i, int j) { return matrixPosition(output,input) * quadraticSize + i * r - i * (i+1) / 2 + j; }

  // buffers for fast evaluation
  double * qiqj;
  double * buffer1;
  void InitBuffers();
  void FreeBuffers();

  int useSingleThread;
  int mkl_max_threads;
  int mkl_dynamic;

  int shallowCopy;
};

#endif

