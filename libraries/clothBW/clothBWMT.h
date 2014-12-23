/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "clothBW" library , Copyright (C) 2014 USC                            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Andy Pierce, Yu Yu Xu, Jernej Barbic                     *
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
 Multi-threaded version of the ClothBW class. 
 It uses the POSIX threads ("pthreads").
*/

#ifndef _CLOTHBWMT_H_
#define _CLOTHBWMT_H_

#include "clothBW.h"

enum ClothBWMT_computationTargetType { STRETCH_SHEAR_FORCE, BEND_FORCE, STRETCH_SHEAR_STIFFNESS, BEND_STIFFNESS };

class ClothBWMT : public ClothBW
{
public:
  
  // constructor that does not require triangleUVs input (it computes UVs automatically; note: the UVs are continuous only within each triangle; the UV map is not global (which is fine, provided one does not want to simulate anisotropic effects) )
  ClothBWMT(int numParticles, double * masses, double * restPositions, int numTriangles,
          int * triangles, int * triangleGroups, int numMaterialGroups,
          double * groupTensileStiffness, double * groupShearStiffness,
          double * groupBendStiffnessU, double * groupBendStiffnessV,
          double * groupDamping, int addGravity=0, int numThreads=1);
  
  // constructor with triangleUVs input
  ClothBWMT(int numParticles, double * masses, double * restPositions, int numTriangles,
          int * triangles, double * traingleUVs, int * triangleGroups, int numMaterialGroups,
          double * groupTensileStiffness, double * groupShearStiffness,
          double * groupBendStiffnessU, double * groupBendStiffnessV,
          double * groupDamping, int addGravity=0, int numThreads=1);	
  
  // constructor from regular clothBW and thread count
  ClothBWMT(ClothBW & clothBW, int numThreads_);
  
  // destructor
  virtual ~ClothBWMT();
  
  // multi-threaded computations
  // compute the internal elastic force, under deformation u
  // note: the force has the same sign as an external force acting on the body (opposite sign as in the StVK class)
  virtual void ComputeForce(double * u, double * f, bool addForce=true); // if addForce is "true", f will be not be reset to zero prior to adding the forces
  
  // compute the damping force
  // unimplemented
  // you can use the damping available in the integrator class
  virtual void ComputeDampingForce(double * u, double * uvel, double * f, bool addForce=false);
  
  virtual void ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix=false);
  
  // for stretch/shear
  int GetStartTriangle(int rank);
  int GetEndTriangle(int rank);
  
  // for bend
  int GetStartQuad(int rank);
  int GetEndQuad(int rank);
  
protected:
  int numThreads;
  
  int * startTriangle;
  int * endTriangle;
  
  int * startQuad;
  int * endQuad;
  
  bool forceAlreadyCleared; // need this to make sure we don't erase a force
  bool matrixAlreadyCleared; // likewise with respect to a stiffness matrix
  
  double * internalForceBuffer;
  SparseMatrix ** sparseMatrixBuffer;
  
  void Initialize();
  void ComputeHelper(enum ClothBWMT_computationTargetType computationTarget, double * u, double * uSecondary, void * target, bool addQuantity);
  
};

#endif

