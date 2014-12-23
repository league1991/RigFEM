/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2014 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
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
#include <math.h>
#include <pthread.h>
#include <vector>
#include <set>
#include "macros.h"
#include "massSpringSystemMT.h"
using namespace std;

//#include "performanceCounter.h"

MassSpringSystemMT::MassSpringSystemMT(int numParticles_, double * masses_, double * restPositions_, int numEdges_, int * edges_, int * edgeGroups_, int numMaterialGroups_, double * groupStiffness_, double * groupDamping_, int addGravity_, int numThreads_) : MassSpringSystem(numParticles_, masses_, restPositions_, numEdges_, edges_, edgeGroups_, numMaterialGroups_, groupStiffness_, groupDamping_, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

MassSpringSystemMT::MassSpringSystemMT(int numParticles_, double * restPositions_, int numQuads, int * quads, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffness, double damping, int addGravity_, int numThreads_): MassSpringSystem(numParticles_, restPositions_, numQuads, quads, surfaceDensity, tensileStiffness, shearStiffness, bendStiffness, damping, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

MassSpringSystemMT::MassSpringSystemMT(int numParticles_, double * restPositions_, MassSpringSystemElementType elementType, int numElements, int * elements, double density, double tensileStiffness, double damping, int addGravity_, int numThreads_): MassSpringSystem(numParticles_, restPositions_, elementType, numElements, elements, density, tensileStiffness, damping, addGravity_), numThreads(numThreads_)
{
  Initialize();
}

MassSpringSystemMT::MassSpringSystemMT(MassSpringSystem & massSpringSystem, int numThreads_): MassSpringSystem(massSpringSystem), numThreads(numThreads_)
{
  Initialize();
}

MassSpringSystemMT::~MassSpringSystemMT()
{
  free(startEdge);
  free(endEdge);
  free(internalForceBuffer);
  for(int i=0; i<numThreads; i++)
    delete(sparseMatrixBuffer[i]);
  free(sparseMatrixBuffer);
}

struct MassSpringSystemMT_threadArg
{
  MassSpringSystemMT * massSpringSystemMT;
  double * u;
  double * uSecondary;
  void * targetBuffer;
  int rank;
  enum MassSpringSystemMT_computationTargetType computationTarget;
};

void * MassSpringSystemMT_WorkerThread(void * arg)
{
  struct MassSpringSystemMT_threadArg * threadArgp = (struct MassSpringSystemMT_threadArg*) arg;
  MassSpringSystemMT * massSpringSystemMT = threadArgp->massSpringSystemMT;
  double * u = threadArgp->u;
  int rank = threadArgp->rank;
  int startEdge = massSpringSystemMT->GetStartEdge(rank);
  int endEdge = massSpringSystemMT->GetEndEdge(rank);

  switch (threadArgp->computationTarget)
  {
    case FORCE: 
    {
      double * targetBuffer = (double*)(threadArgp->targetBuffer);
      massSpringSystemMT->AddForce(u, targetBuffer, startEdge, endEdge);
    }
    break;

    case DAMPINGFORCE:
    {
      double * uvel = u;
      double * targetBuffer = (double*)(threadArgp->targetBuffer);
      massSpringSystemMT->AddDampingForce(uvel, targetBuffer, startEdge, endEdge);
    }
    break;

    case STIFFNESSMATRIX:
    {
      SparseMatrix * targetBuffer = (SparseMatrix*)(threadArgp->targetBuffer);
      massSpringSystemMT->AddStiffnessMatrix(u, targetBuffer, startEdge, endEdge);
    }
    break;

    case HESSIANAPPROXIMATION:
    {
      SparseMatrix * targetBuffer = (SparseMatrix*)(threadArgp->targetBuffer);
      double * uSecondary = threadArgp->uSecondary;
      massSpringSystemMT->AddHessianApproximation(u, uSecondary, targetBuffer, startEdge, endEdge);
    }
    break;

    default:
      printf("Error: unknown computation type.\n");
      exit(1);
    break;
  }

  return NULL;
}

void MassSpringSystemMT::Initialize()
{
  internalForceBuffer = (double*) malloc (sizeof(double) * numThreads * 3 * numParticles);

  // generate skeleton matrices
  sparseMatrixBuffer = (SparseMatrix**) malloc (sizeof(SparseMatrix*) * numThreads);

  SparseMatrix * sparseMatrix;
  GetStiffnessMatrixTopology(&sparseMatrix);
  for(int i=0; i<numThreads; i++)
    sparseMatrixBuffer[i] = new SparseMatrix(*sparseMatrix);

  // split the workload
  startEdge = (int*) malloc (sizeof(int) * numThreads);
  endEdge = (int*) malloc (sizeof(int) * numThreads);

  int remainder = numEdges % numThreads;
  // the first 'remainder' nodes will process one edge more
  int jobSize = numEdges / numThreads;

  for(int rank=0; rank < numThreads; rank++)
  {
    if (rank < remainder)
    {
      startEdge[rank] = rank * (jobSize+1);
      endEdge[rank] = (rank+1) * (jobSize+1);
    }
    else
    {
      startEdge[rank] = remainder * (jobSize+1) + (rank-remainder) * jobSize;
      endEdge[rank] = remainder * (jobSize+1) + ((rank-remainder)+1) * jobSize;
    }
  }

  delete(sparseMatrix);
  printf("Total edges: %d \n",numEdges);
  printf("Num threads: %d \n",numThreads);
  printf("Canonical job size: %d \n",jobSize);
  printf("Num threads with job size augmented by one edge: %d \n",remainder);
}

void MassSpringSystemMT::ComputeHelper(enum MassSpringSystemMT_computationTargetType computationTarget, double * u, double * uSecondary, void * target, bool addQuantity)
{
  // launch threads
  int numParticles3 = 3*numParticles;
  struct MassSpringSystemMT_threadArg * threadArgv = (struct MassSpringSystemMT_threadArg*) malloc (sizeof(struct MassSpringSystemMT_threadArg) * numThreads);

  pthread_t * tid = (pthread_t*) malloc (sizeof(pthread_t) * numThreads);

  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].massSpringSystemMT = this;
    threadArgv[i].u = u;
    threadArgv[i].rank = i;
    threadArgv[i].computationTarget = computationTarget;
  }

  switch(computationTarget)
  {
    case FORCE:
    case DAMPINGFORCE:
    {
      for(int i=0; i<numThreads; i++)
        threadArgv[i].targetBuffer = (void*)(&internalForceBuffer[i * numParticles3]);
      memset(internalForceBuffer, 0, sizeof(double) * numParticles3 * numThreads);
    }
    break;

    case STIFFNESSMATRIX:
    case HESSIANAPPROXIMATION:
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
    if (pthread_create(&tid[i], NULL, MassSpringSystemMT_WorkerThread, &threadArgv[i]) != 0)
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
    case FORCE:
    case DAMPINGFORCE:
    {
      double * f = (double*) target;

      if (!addQuantity)
        memset(f, 0, sizeof(double) * numParticles3);

      for(int i=0; i<numThreads; i++)
      {
        double * source = &internalForceBuffer[i * numParticles3];
        for(int j=0; j<numParticles3; j++)
          f[j] += source[j];
      }

      if ((computationTarget == FORCE) && addGravity)
        ComputeGravity(f, true);
    }
    break;

    case STIFFNESSMATRIX:
    case HESSIANAPPROXIMATION:
    {
      SparseMatrix * targetK = (SparseMatrix*) target;

      if (!addQuantity)
        targetK->ResetToZero();

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

void MassSpringSystemMT::ComputeForce(double * u, double * f, bool addForce)
{
  //PerformanceCounter forceCounter;
  ComputeHelper(FORCE, u, NULL, (void*)f, addForce);

  //forceCounter.StopCounter();
  //printf("Internal forces: %G\n", forceCounter.GetElapsedTime());
}

void MassSpringSystemMT::ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix)
{
  //PerformanceCounter stiffnessCounter;
  ComputeHelper(STIFFNESSMATRIX, u, NULL, (void*)K , addMatrix);

  //stiffnessCounter.StopCounter();
  //printf("Stiffness matrix: %G\n", stiffnessCounter.GetElapsedTime());
}

void MassSpringSystemMT::ComputeDampingForce(double * uvel, double * f, bool addForce)
{
  ComputeHelper(DAMPINGFORCE, uvel, NULL, (void*)f, addForce);
}

void MassSpringSystemMT::ComputeHessianApproximation(double * u, double * du, SparseMatrix * dK, bool addMatrix)
{
  ComputeHelper(HESSIANAPPROXIMATION, u, du, (void*) dK, addMatrix);
}

int MassSpringSystemMT::GetStartEdge(int rank)
{
  return startEdge[rank];
}

int MassSpringSystemMT::GetEndEdge(int rank)
{
  return endEdge[rank];
}

