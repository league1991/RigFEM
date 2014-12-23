/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2014 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic, Daniel Schroeder                          *
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

#ifndef _MASS_SPRING_SYSTEM_MT_H_
#define _MASS_SPRING_SYSTEM_MT_H_

#include "massSpringSystem.h"

/*
   Multi-threaded version of the MassSpringSystem class. 
   It uses the POSIX threads ("pthreads").
   See also massSpringSystem.h
*/

enum MassSpringSystemMT_computationTargetType { FORCE, STIFFNESSMATRIX, DAMPINGFORCE, HESSIANAPPROXIMATION };

class MassSpringSystemMT : public MassSpringSystem
{
public:

  // creates a mass spring from scratch
  MassSpringSystemMT(int numParticles, double * masses, double * restPositions,
    int numEdges, int * edges, int * edgeGroups, int numMaterialGroups, double * groupStiffness, double * groupDamping, int addGravity=0, int numThreads=1);

  MassSpringSystemMT(int numParticles, double * restPositions, int numQuads, int * quads, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffness, double damping, int addGravity=0, int numThreads=1); // creates the mass spring system from a quad surface mesh (cloth)

  MassSpringSystemMT(int numParticles, double * restPositions, MassSpringSystemElementType elementType, int numElements, int * elements, double density, double tensileStiffness, double damping, int addGravity=0, int numThreads=1); // creates the mass spring system from a tet meh

  MassSpringSystemMT(MassSpringSystem & massSpringSystem, int numThreads);

  virtual ~MassSpringSystemMT();

  virtual void ComputeForce(double * u, double * f, bool addForce=false); 
  virtual void ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix=false);
  virtual void ComputeDampingForce(double * uvel, double * f, bool addForce=false); 
  // computes an approximation to dK, assuming the deformations change from u to u + du
  virtual void ComputeHessianApproximation(double * u, double * du, SparseMatrix * dK, bool addMatrix=false);

  int GetStartEdge(int rank);
  int GetEndEdge(int rank);

protected:
  int numThreads;
  int * startEdge, * endEdge;
  double * internalForceBuffer;
  SparseMatrix ** sparseMatrixBuffer;

  void Initialize();
  void ComputeHelper(enum MassSpringSystemMT_computationTargetType computationTarget, double * u, double * uSecondary, void * target, bool addQuantity);

};

#endif

