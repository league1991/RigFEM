/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "isotropic hyperelastic FEM" library , Copyright (C) 2014 USC         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin                            *
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
#include "isotropicHyperelasticFEMMT.h"

IsotropicHyperelasticFEMMT::IsotropicHyperelasticFEMMT(TetMesh * tetMesh_, IsotropicMaterial * isotropicMaterial_, double principalStretchThreshold_, bool addGravity_, double g_, int numThreads_) :
  IsotropicHyperelasticFEM(tetMesh_, isotropicMaterial_, principalStretchThreshold_, addGravity_, g_),
  numThreads(numThreads_)
{
  Initialize();
}

IsotropicHyperelasticFEMMT::~IsotropicHyperelasticFEMMT()
{
  free(startElement);
  free(endElement);
  free(energyBuffer);
  free(internalForceBuffer);
  for(int i=0; i<numThreads; i++)
    delete(tangentStiffnessMatrixBuffer[i]);
  free(tangentStiffnessMatrixBuffer);
}

struct IsotropicHyperelasticFEMMT_threadArg
{
  IsotropicHyperelasticFEMMT * isotropicHyperelasticFEMMT;
  double * u;
  double * energy;
  double * internalForces;
  SparseMatrix *  tangentStiffnessMatrix;
  int computationMode;
  int rank;
  int exitCode;
};

void * IsotropicHyperelasticFEMMT_WorkerThread(void * arg)
{
  struct IsotropicHyperelasticFEMMT_threadArg * threadArgp = (struct IsotropicHyperelasticFEMMT_threadArg*) arg;
  IsotropicHyperelasticFEMMT * isotropicHyperelasticFEMMT = threadArgp->isotropicHyperelasticFEMMT;
  double * u = threadArgp->u;
  double * energy = threadArgp->energy;
  double * internalForces = threadArgp->internalForces;
  SparseMatrix * tangentStiffnessMatrix = threadArgp->tangentStiffnessMatrix;
  int computationMode = threadArgp->computationMode;
  int rank = threadArgp->rank;
  int startElement = isotropicHyperelasticFEMMT->GetStartElement(rank);
  int endElement = isotropicHyperelasticFEMMT->GetEndElement(rank);

  //printf("%d %d\n", startElement, endElement);
  int code = isotropicHyperelasticFEMMT->GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(startElement, endElement, u, energy, internalForces, tangentStiffnessMatrix, computationMode); 
  threadArgp->exitCode = code;

  return NULL;
}

void IsotropicHyperelasticFEMMT::Initialize()
{
  energyBuffer = (double*) malloc (sizeof(double) * numThreads);
  internalForceBuffer = (double*) malloc (sizeof(double) * numThreads * 3 * tetMesh->getNumVertices());

  // generate skeleton matrices
  tangentStiffnessMatrixBuffer = (SparseMatrix**) malloc (sizeof(SparseMatrix*) * numThreads);

  SparseMatrix * sparseMatrix;
  GetStiffnessMatrixTopology(&sparseMatrix);
  for(int i=0; i<numThreads; i++)
    tangentStiffnessMatrixBuffer[i] = new SparseMatrix(*sparseMatrix);

  // split the workload
  int numElements = tetMesh->getNumElements();
  startElement = (int*) malloc (sizeof(int) * numThreads);
  endElement = (int*) malloc (sizeof(int) * numThreads);

  int remainder = numElements % numThreads;
  // the first 'remainder' nodes will process one edge more
  int jobSize = numElements / numThreads;

  for(int rank=0; rank < numThreads; rank++)
  {
    if (rank < remainder)
    {
      startElement[rank] = rank * (jobSize+1);
      endElement[rank] = (rank+1) * (jobSize+1);
    }
    else
    {
      startElement[rank] = remainder * (jobSize+1) + (rank-remainder) * jobSize;
      endElement[rank] = remainder * (jobSize+1) + ((rank-remainder)+1) * jobSize;
    }
  }

  printf("Total elements: %d \n", numElements);
  printf("Num threads: %d \n", numThreads);
  printf("Canonical job size: %d \n", jobSize);
  printf("Num threads with job size augmented by one edge: %d \n", remainder);
}

int IsotropicHyperelasticFEMMT::GetEnergyAndForceAndTangentStiffnessMatrixHelper(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode)
{
  GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(u, energy, internalForces, tangentStiffnessMatrix, computationMode);

  // launch threads
  struct IsotropicHyperelasticFEMMT_threadArg * threadArgv = (struct IsotropicHyperelasticFEMMT_threadArg*) malloc (sizeof(struct IsotropicHyperelasticFEMMT_threadArg) * numThreads);

  pthread_t * tid = (pthread_t*) malloc (sizeof(pthread_t) * numThreads);

  int numVertices3 = 3 * tetMesh->getNumVertices();
  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].isotropicHyperelasticFEMMT = this;
    threadArgv[i].u = u;
    threadArgv[i].energy = &energyBuffer[i];
    threadArgv[i].internalForces = &internalForceBuffer[i * numVertices3];
    threadArgv[i].tangentStiffnessMatrix = tangentStiffnessMatrixBuffer[i];
    threadArgv[i].computationMode = computationMode;
    threadArgv[i].rank = i;
  }

  // clear internal buffers
  memset(energyBuffer, 0, sizeof(double) * numThreads);
  memset(internalForceBuffer, 0, sizeof(double) * numThreads * numVertices3);
  for(int i=0; i<numThreads; i++)  
    tangentStiffnessMatrixBuffer[i]->ResetToZero();

  for(int i=0; i<numThreads; i++)  
  {
    if (pthread_create(&tid[i], NULL, IsotropicHyperelasticFEMMT_WorkerThread, &threadArgv[i]) != 0)
    {
      printf("Error: unable to launch thread %d.\n", i);
      return 1;
    }
  }

  int code = 0;
  for(int i=0; i<numThreads; i++)
  {
    if (pthread_join(tid[i], NULL) != 0)
    {
      printf("Error: unable to join thread %d.\n", i);
      exit(1);
    }
    if (threadArgv[i].exitCode != 0)
      code = 1;
  }

  free(threadArgv);
  free(tid);

  for(int i=0; i<numThreads; i++)
  {
    if (computationMode & COMPUTE_ENERGY)
      *energy += energyBuffer[i];

    if (computationMode & COMPUTE_INTERNALFORCES)
    {
      double * source = &internalForceBuffer[i * numVertices3];
      for(int j=0; j<numVertices3; j++)
        internalForces[j] += source[j];
    }

    if (computationMode & COMPUTE_TANGENTSTIFFNESSMATRIX)
    {
      *tangentStiffnessMatrix += *(tangentStiffnessMatrixBuffer[i]);
    }
  }

  return code;
}

int IsotropicHyperelasticFEMMT::GetStartElement(int rank)
{
  return startElement[rank];
}

int IsotropicHyperelasticFEMMT::GetEndElement(int rank)
{
  return endElement[rank];
}

