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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clothBWMT.h"

#include <pthread.h>

ClothBWMT::ClothBWMT(int numParticles_, double * masses_, double * restPositions_, int numTriangles_,
                      int * triangles_, int * triangleGroups_, int numMaterialGroups_,
                       double * groupTensileStiffness_, double * groupShearStiffness_,
                       double * groupBendStiffnessU_, double * groupBendStiffnessV_,
                       double * groupDamping_, int addGravity_, int numThreads_) : 
                            ClothBW(numParticles_, masses_, restPositions_, numTriangles_,
                                    triangles_, triangleGroups_, numMaterialGroups_,
                                    groupTensileStiffness_, groupShearStiffness_,
                                    groupBendStiffnessU_, groupBendStiffnessV_,
                                    groupDamping_, addGravity_), numThreads(numThreads_)
                                                                                        
{
  Initialize();
}

// constructor with triangleUVs input
ClothBWMT::ClothBWMT(int numParticles_, double * masses_, double * restPositions_, int numTriangles_,
           int * triangles_, double * triangleUVs_, int * triangleGroups_, int numMaterialGroups_,
           double * groupTensileStiffness_, double * groupShearStiffness_,
           double * groupBendStiffnessU_, double * groupBendStiffnessV_,
           double * groupDamping_, int addGravity_, int numThreads_) : 
                          ClothBW(numParticles_, masses_, restPositions_, numTriangles_,
                                  triangles_, triangleUVs_, triangleGroups_, numMaterialGroups_,
                                  groupTensileStiffness_, groupShearStiffness_,
                                  groupBendStiffnessU_, groupBendStiffnessV_,
                                  groupDamping_, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

ClothBWMT::ClothBWMT(ClothBW & clothBW, int numThreads_) : ClothBW(clothBW), numThreads(numThreads_)
{
  Initialize();
}

// destructor
ClothBWMT::~ClothBWMT()
{
  free(startTriangle);
  free(endTriangle);
  
  free(startQuad);
  free(endQuad);
  
  free(internalForceBuffer);
}

// data structure to hold data for each thread
struct ClothBWMT_threadArg
{
  ClothBWMT * clothBW_MT;
  double * u;
  double * uSecondary;
  void * targetBuffer;
  int rank;
  enum ClothBWMT_computationTargetType computationTarget;
};

// global function
void * ClothBWMT_WorkerThread(void * arg)
{
  // cast to struct
  struct ClothBWMT_threadArg * threadArgp = (struct ClothBWMT_threadArg*) arg;
  
  // copy info into local vars
  ClothBWMT * clothBW_MT = threadArgp->clothBW_MT;
  double * u = threadArgp->u;
  int rank = threadArgp->rank;
  int startTriangle = clothBW_MT->GetStartTriangle(rank);
  int endTriangle = clothBW_MT->GetEndTriangle(rank);
  int startQuad = clothBW_MT->GetStartQuad(rank);
  int endQuad = clothBW_MT->GetEndQuad(rank);
  
  // now put the thread to work
  switch (threadArgp->computationTarget)
  {
    case STRETCH_SHEAR_FORCE:
    {
      double * targetBuffer = (double*)(threadArgp->targetBuffer);
      clothBW_MT->AddStretchAndShearForce(u, targetBuffer, startTriangle, endTriangle);
    }
      break;
      
    case BEND_FORCE:
    {
      double * targetBuffer = (double*)(threadArgp->targetBuffer);
      clothBW_MT->AddBendForce(u, targetBuffer, startQuad, endQuad);
    }
      break;
      
    case STRETCH_SHEAR_STIFFNESS:
    {
      SparseMatrix * targetBuffer = (SparseMatrix*)(threadArgp->targetBuffer);
      clothBW_MT->AddStretchAndShearStiffnessMatrix(u, targetBuffer, startTriangle, endTriangle);
    }
      break;
      
    case BEND_STIFFNESS:
    {
      SparseMatrix * targetBuffer = (SparseMatrix*)(threadArgp->targetBuffer);
      clothBW_MT->AddBendStiffnessMatrix(u, targetBuffer, startQuad, endQuad);
    }
      break;
      
    default:
      printf("Error: unknown computation type.\n");
      exit(1);
      break;
  }
  
  return NULL;
  
}

// multi-threaded computations
// compute the internal elastic force, under deformation u
// note: the force has the same sign as an external force acting on the body (opposite sign as in the StVK class)
// if addForce is "true", f will be not be reset to zero prior to adding the forces
void ClothBWMT::ComputeForce(double * u, double * f, bool addForce)
{
  forceAlreadyCleared = false;
  
  if (_computationConditions[0])
    ComputeHelper(STRETCH_SHEAR_FORCE, u, NULL, (void*)f, addForce);
  if (_computationConditions[1])
    ComputeHelper(BEND_FORCE, u, NULL, (void*)f, addForce);
}

// compute the damping force
// unimplemented
void ClothBWMT::ComputeDampingForce(double * u, double * uvel, double * f, bool addForce)
{
  
}

void ClothBWMT::ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix)
{
  //ComputeHelper(STIFFNESSMATRIX, u, NULL, (void*)K, addMatrix);
  matrixAlreadyCleared = false;
  
  if (_computationConditions[2])
    ComputeHelper(STRETCH_SHEAR_STIFFNESS, u, NULL, (void*)K, addMatrix);
  
  if (_computationConditions[3])
    ComputeHelper(BEND_STIFFNESS, u, NULL, (void*)K, addMatrix);
}

// for stretch/shear
int ClothBWMT::GetStartTriangle(int rank)
{
  return startTriangle[rank];
}

int ClothBWMT::GetEndTriangle(int rank)
{
  return endTriangle[rank];
}

// for bend
int ClothBWMT::GetStartQuad(int rank)
{
  return startQuad[rank];
}

int ClothBWMT::GetEndQuad(int rank)
{
  return endQuad[rank];
}

void ClothBWMT::Initialize()
{
  // size of our force buffer (3 x numParticles for each thread)
  internalForceBuffer = (double*) malloc (sizeof(double) * numThreads * 3 * numParticles);
  
  // generate skeleton matrices
  sparseMatrixBuffer = (SparseMatrix**) malloc (sizeof(SparseMatrix*) * numThreads);
  
  SparseMatrix * sparseMatrix;
  GenerateStiffnessMatrixTopology(&sparseMatrix);
  for(int i=0; i<numThreads; i++)
    sparseMatrixBuffer[i] = new SparseMatrix(*sparseMatrix);
  
  // split the workload
  startTriangle = (int*) malloc (sizeof(int) * numThreads);
  endTriangle = (int*) malloc (sizeof(int) * numThreads);
  
  startQuad = (int*) malloc (sizeof(int) * numThreads);
  endQuad = (int*) malloc (sizeof(int) * numThreads);
  
  int remainderTriangles = numTriangles % numThreads;
  int remainderQuads = numQuads % numThreads;
  
  // the first 'remainder' nodes will process one triangle + one quad more
  
  int jobSizeTriangles = numTriangles / numThreads;
  int jobSizeQuads = numQuads / numThreads;
  
  for(int rank=0; rank < numThreads; rank++)
  {
    // triangles
    if (rank < remainderTriangles)
    {
      startTriangle[rank] = rank * (jobSizeTriangles+1);
      endTriangle[rank] = (rank+1) * (jobSizeTriangles+1);
    }
    else
    {
      startTriangle[rank] = remainderTriangles * (jobSizeTriangles+1) + (rank-remainderTriangles) * jobSizeTriangles;
      endTriangle[rank] = remainderTriangles * (jobSizeTriangles+1) + ((rank-remainderTriangles)+1) * jobSizeTriangles;
    }
    
    // quads
    if (rank < remainderQuads)
    {
      startQuad[rank] = rank * (jobSizeQuads+1);
      endQuad[rank] = (rank+1) * (jobSizeQuads+1);
    }
    else
    {
      startQuad[rank] = remainderQuads * (jobSizeQuads+1) + (rank-remainderQuads) * jobSizeQuads;
      endQuad[rank] = remainderQuads * (jobSizeQuads+1) + ((rank-remainderQuads)+1) * jobSizeQuads;
    }
    
  }
  
  printf("Total triangles: %d \n",numTriangles);
  printf("Total quads: %d \n",numQuads);
  
  printf("Num threads: %d \n",numThreads);
  printf("Canonical job size (triangles): %d \n",jobSizeTriangles);
  printf("Canonical job size (quads): %d \n",jobSizeQuads);
  
  printf("Num threads with job size augmented by one triangle: %d \n",remainderTriangles);
  printf("Num threads with job size augmented by one quad: %d \n",remainderQuads);
  
}

void ClothBWMT::ComputeHelper(enum ClothBWMT_computationTargetType computationTarget, double * u, double * uSecondary, void * target, bool addQuantity)
{
  // launch threads
  int numParticles3 = 3*numParticles;
  struct ClothBWMT_threadArg * threadArgv = (struct ClothBWMT_threadArg*) malloc (sizeof(struct ClothBWMT_threadArg) * numThreads);
  
  pthread_t * tid = (pthread_t*) malloc (sizeof(pthread_t) * numThreads);
  
  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].clothBW_MT = this;
    threadArgv[i].u = u;
    threadArgv[i].rank = i;
    threadArgv[i].computationTarget = computationTarget;
  }
  
  switch(computationTarget)
  {
    // all force calcs need same size buffer, so we just fall through here
    case STRETCH_SHEAR_FORCE:
    case BEND_FORCE:
    {
      for(int i=0; i<numThreads; i++)
        threadArgv[i].targetBuffer = (void*)(&internalForceBuffer[i * numParticles3]);
      memset(internalForceBuffer, 0, sizeof(double) * numParticles3 * numThreads);
    }
      break;
      
    // likewise for the stiffness matrix calculations
    case STRETCH_SHEAR_STIFFNESS:
    case BEND_STIFFNESS:
    {
      for(int i=0; i<numThreads; i++)
      {
        threadArgv[i].targetBuffer = (void*)(sparseMatrixBuffer[i]);
        sparseMatrixBuffer[i]->ResetToZero();
      }
    }
      break;
      
    default:
      printf("Error: unknown computation type.\n");
      exit(1);
      break;
  }
  
  for(int i=0; i<numThreads; i++)  
  {
    if (pthread_create(&tid[i], NULL, ClothBWMT_WorkerThread, &threadArgv[i]) != 0)
    {
      printf("Error: unable to launch thread %d.\n", i);
      exit(1);
    }
  }
  
  for(int i=0; i<numThreads; i++)
  {
    if (pthread_join(tid[i], NULL) != 0)
    {
      printf("Error: unable to join thread %d.\n", i);
      exit(1);
    }
  }
  
  free(threadArgv);
  free(tid);
  
  // assemble results
  switch(computationTarget)
  {
    case STRETCH_SHEAR_FORCE:   
    case BEND_FORCE:
    {
      double * f = (double*) target;
      
      if (!addQuantity && !forceAlreadyCleared)
      {
        memset(f, 0, sizeof(double) * numParticles3);
        forceAlreadyCleared = true;
      }
      
      for(int i=0; i<numThreads; i++)
      {
        double * source = &internalForceBuffer[i * numParticles3];
        for(int j=0; j<numParticles3; j++)
          f[j] += source[j];
      }
      
      if (addGravity)
        ComputeGravity(f, true);
    }
      break;
      
    case STRETCH_SHEAR_STIFFNESS:
    case BEND_STIFFNESS:
    {
      SparseMatrix * targetK = (SparseMatrix*) target;
      
      if (!addQuantity && !matrixAlreadyCleared)
      {
        targetK->ResetToZero();
        matrixAlreadyCleared = true;
      }
      
      for(int i=0; i<numThreads; i++)
      {
        *targetK += *(sparseMatrixBuffer[i]);
      }
    }
      break;
      
    default:
      printf("Error: unknown computation type.\n");
      exit(1);
      break;
  }
}

