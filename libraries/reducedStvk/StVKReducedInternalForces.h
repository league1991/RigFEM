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
  This class computes the reduced internal forces for a reduced St. Venant Kirchhoff deformable model.
  Each element of the force vector is a cubic polynomial in the reduced coordinates q.
  The force vector and reduced coordinates are vectors of length r, 
  where r is the dimension of the reduced basis.

  The class can perform the following:
  1. Precompute the coefficients of the cubic polynomials
  2. Given the coefficients, evaluate the reduced internal forces quickly. 
     The evaluation routine uses BLAS Level 2 and 3 routines. 
     The evaluation routine has been optimized extensively, and is very efficient.

  See also StVKInternalForces.h .
*/

#ifndef _STVKREDUCEDINTERNALFORCES_H_
#define _STVKREDUCEDINTERNALFORCES_H_

#include "StVKInternalForces.h"

class StVKReducedInternalForces
{
public:

  // compute the reduced cubic coefficients from scratch, for the given basis and volumetric mesh
  // r : dimension of the basis
  // U : the basis matrix (it must be of size 3*n x r, where n is the number of vertices in the volumetric mesh)
  // volumetricMesh, precomputedABCDIntegrals: same meaning as in StVKInternalForces.h
  // initOnly:
  //   0: normal mode (default, recommended)
  //   1: do not actually run the computation (advanced low-level mode; must then manually call ProcessElements)
  // note: this computation can last for several minutes, depending on r and n
  StVKReducedInternalForces(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, int initOnly=0, bool addGravity=false, double g=9.81, int verbose=1);

  // load the previously computed coefficients from a file 
  // if r>=0 is specified, only up to first r modes will be used; if r=-1 (default), all modes will be used
  // bigEndianMachine allows you to load little endian data (e.g. PC, Mac Intel) on a big endian machine (e.g. Mac PowerPC); default: 0 (no conversion)
  StVKReducedInternalForces(const char * filename, int r=-1, int bigEndianMachine=0, int verbose=1);
  StVKReducedInternalForces(FILE * fin, int r=-1, int bigEndianMachine=0, int verbose=1); // read from binary stream

  ~StVKReducedInternalForces();

  // saves coefficients to a disk file (in binary format, the convention is to use the .cub file extension)
  int Save(const char * filename);
  int Save(FILE * fout); // saves to stream
  
  // evaluates the reduced internal forces for the given configuration q, result is written into fq (which must be a pre-allocated vector of length r)
  // see also the f_int(x) comment in StVKInternalForces.h on the sign of fq; the comment applies here too
  void Evaluate(double * q, double * fq);
  double EvaluateComponent(double * q, int componentIndex); // evaluates just one force component (a scalar)
  void EvaluateLinear(double * q, double * fq); // evaluates just the linear terms (constant terms are zero in all polynomials as reduced internal force is zero at the origin)

  // enables or disables the gravity (note: you can also set this in the constructor; use this routine to turn the gravity on/off during the simulation)
  void SetGravity(bool addGravity, double g, VolumetricMesh * volumetricMesh_=NULL, double * U_=NULL) { this->addGravity = addGravity; this->g=g; InitGravity(volumetricMesh_, U_); } // if AddGravity is enabled, Evaluate will subtract the gravity force from the reduced internal forces (note: subtraction, not addition, is used; see the comment in StVKInternalForces.h, which also applies here) 

  inline int Getr() { return r; }
  // report r from the given polynomial binary file 
  static int GetrFromFile(const char * filename);

  void Scale(double scalingFactor); // scales all coefficients with the given scaling factor; i.e., linearly and uniformly scale the stiffness of the model; frequency spectrum will scale linearly by sqrt(scalingFactor)

  // ==== all routines below here are low-level research routines (expected to be rarely used) ===

  // explicitly query coefficients of polynomials
  // "index" denotes the index of the queried polynomial; 0 <= index < r
  inline double linearCoef(int index, int i) { return linearCoef_[linearCoefPos(index,i)]; } // coefficient at q_i; must pass 0 <= i < r
  inline double quadraticCoef(int index, int i, int j) { return quadraticCoef_[quadraticCoefPos(index,i,j)]; } // coefficient at q_i q_j; must pass 0 <= i <= j < r
  inline double cubicCoef(int index, int i,int j, int k) { return cubicCoef_[cubicCoefPos(index,i,j,k)]; } // coefficient at q_i q_j q_k; must pass 0 <= i <= j <= k < r

  // prints the polynomial to standard output
  void PrintPolynomial(); // in Mathematica format
  void PrintLinearCoefficients(); // raw data; one coefficient per line
  void PrintQuadraticCoefficients(); // raw data; one coefficient per line
  void PrintCubicCoefficients(); // raw data; one coefficient per line

  // compares the value of the polynomial with evaluation via U^T f(Uq)
  // prints out the two values (they should be equal)
  void TestPolynomial(double * q, StVKInternalForces * stVKInternalForces, const char * modalMatrixFilename);

  // sorts three integers in ascending order
  static void tripleSort(int & a, int & b, int & c);

  // build histograms (for statistics)
  // x-scale is log10, will fill buckets from 10^{lowExp} to 10^{highExp}
  // warning will be generated if elements fall out of the given range
  // histogram must be pre-allocated array
  int CubicCoefficientsHistogramLog10(int * histogram,int lowExp, int highExp);

  // finds max coefficient in terms of abs value among all terms
  double MaxAbsCoefficient();
  double MaxAbsLinearCoefficient();
  double MaxAbsQuadraticCoefficient();
  double MaxAbsCubicCoefficient();
  double AverageAbsCubicCoefficient();

  static inline int GetLinearSize(int r) { return r; } // number of double-precision entries
  static inline int GetQuadraticSize(int r) { return r*(r+1)/2; } // number of double-precision entries
  static inline int GetCubicSize(int r) { return r * (r+1) * (r+2) / 6; } // number of double-precision entries
  static inline int GetTotalFileSize(int r) { return r * (GetLinearSize(r) + GetQuadraticSize(r) + GetCubicSize(r)) * sizeof(double) + 4 * sizeof(int); } // in bytes

  inline int GetLinearSize() { return GetLinearSize(r); }
  inline int GetQuadraticSize() { return GetQuadraticSize(r); }
  inline int GetCubicSize() { return GetCubicSize(r); }

  // returns pointers to internal storage of coefficients
  inline double * GetLinearTermsBuffer() { return linearCoef_; }
  inline double * GetQuadraticTermsBuffer() { return quadraticCoef_; }
  inline double * GetCubicTermsBuffer() { return cubicCoef_; }

  // computes the contributions from voxels startElements to endElement-1
  void ProcessElements(int startElement, int endElement, double ** target=NULL);

  void UseSingleThread(int useSingleThread);

  // makes shallow copies of all pointers, except those initialized by InitBuffers
  // use this if you want to Evaluate two or more identical models (i.e., two copies of an object) in parallel (to ensure thread safety)
  // you do not need to use this if you are Evaluating a single model in parallel (e.g., using the MT derived class)
  StVKReducedInternalForces * ShallowClone();

  // saves a model with r=0
  static int SaveEmptyCub(const char * filename);
  static int SaveEmptyCub(FILE * fout);

protected:

  // Layout of polynomial coefficients in memory:

  // linear terms:
  // r
  // r
  // r
  // r
  // r
  // r
  // r

  //Each "r" is a vector of length r.

  // quadratic terms:
  // q
  // q
  // q
  // q
  // q
  // q
  // q
  
  // Each "q" is an upper triangle of vectors of length r:
  // rrrrrrr
  //  rrrrrr
  //   rrrrr
  //    rrrr
  //     rrr
  //      rr
  //       r

  // cubic terms:
  // c
  // c
  // c
  // c
  // c
  // c
  // c
  
  // Each "c" is a 3D upper triangular pyramid of coefficients
  // (r*(r+1)*(r+2) / 6 coefficients total)

  VolumetricMesh * volumetricMesh;
  StVKElementABCD * precomputedIntegrals;
  int numElementVertices;
  double * lambdaLame, * muLame;

  double * unitReducedGravityForce;
  double * reducedGravityForce;
  bool addGravity;
  double g;
  void InitGravity(VolumetricMesh * volumetricMesh_=NULL, double * U_=NULL);

  void InitComputation(int r, double * U, VolumetricMesh * volumetricMesh);

  void GetSizes(int r, int * linearSize, int * quadraticSize, int * cubicSize);
  void SetSizes(int r);

  int n,r,r2;
  double * U;

  double * linearCoef_;
  double * quadraticCoef_;
  double * cubicCoef_;

  int linearSize;
  int quadraticSize;
  int cubicSize;

  // both index and i are free to be anywhere on 0...r-1
  inline int linearCoefPos(int index, int i) {return index * linearSize + i;}
  // assumes i<=j, index is free
  inline int quadraticCoefPos(int index, int i, int j) { return index * quadraticSize + i * r - i * (i+1) / 2 + j; }
  // assumes i<=j<=k, index is free
  inline int cubicCoefPos(int index, int i, int j, int k);

  void little2big(void * input, void * output, int numBytes);

  // acceleration data
  //double * qij;
  double * qiqj;
  void InitBuffers();
  void FreeBuffers();

  int useSingleThread;
  int mkl_max_threads;
  int mkl_dynamic;
  int shallowCopy;

  int verbose;

  int LoadFromStream(FILE * fin, int rTarget, int bigEndianMachine);
};

inline int StVKReducedInternalForces::cubicCoefPos(int index, int i, int j, int k)
{
  return index * cubicSize +
  + (-j + (i-j)*(-1+i+j)/2 + k + (j-i)*r + i * (2 - 3*i + i*i + 6 * r - 3 * i * r + 3 * r2 ) / 6);
}

#endif

