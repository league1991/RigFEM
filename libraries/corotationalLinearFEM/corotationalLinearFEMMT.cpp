/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "corotational linear FEM" library , Copyright (C) 2014 USC            *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <vector>
#include <set>
#include "macros.h"
#include "corotationalLinearFEMMT.h"
using namespace std;

CorotationalLinearFEMMT::CorotationalLinearFEMMT(TetMesh * tetMesh, int numThreads_) : CorotationalLinearFEM(tetMesh), numThreads(numThreads_)
{
  Initialize();
}

CorotationalLinearFEMMT::~CorotationalLinearFEMMT()
{
  free(startElement);
  free(endElement);
  free(internalForceBuffer);
  for(int i=0; i<numThreads; i++)
    delete(stiffnessMatrixBuffer[i]);
  free(stiffnessMatrixBuffer);
}

struct CorotationalLinearFEMMT_threadArg
{
  CorotationalLinearFEMMT * corotationalLinearFEMMT;
  double * u;
  double * f;
  SparseMatrix * stiffnessMatrix;
  int warp;
  int rank;
};

void * CorotationalLinearFEMMT_WorkerThread(void * arg)
{
  struct CorotationalLinearFEMMT_threadArg * threadArgp = (struct CorotationalLinearFEMMT_threadArg*) arg;
  CorotationalLinearFEMMT * corotationalLinearFEMMT = threadArgp->corotationalLinearFEMMT;
  double * u = threadArgp->u;
  double * f = threadArgp->f;
  SparseMatrix * stiffnessMatrix = threadArgp->stiffnessMatrix;
  int warp = threadArgp->warp;
  int rank = threadArgp->rank;
  int startElement = corotationalLinearFEMMT->GetStartElement(rank);
  int endElement = corotationalLinearFEMMT->GetEndElement(rank);

  //printf("%d %d\n", startElement, endElement);
  corotationalLinearFEMMT->ComputeForceAndStiffnessMatrixOfSubmesh(u, f, stiffnessMatrix, warp, startElement, endElement);

  return NULL;
}

void CorotationalLinearFEMMT::Initialize()
{
  internalForceBuffer = (double*) malloc (sizeof(double) * numThreads * 3 * tetMesh->getNumVertices());

  // generate skeleton matrices
  stiffnessMatrixBuffer = (SparseMatrix**) malloc (sizeof(SparseMatrix*) * numThreads);

  SparseMatrix * sparseMatrix;
  GetStiffnessMatrixTopology(&sparseMatrix);
  for(int i=0; i<numThreads; i++)
    stiffnessMatrixBuffer[i] = new SparseMatrix(*sparseMatrix);

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

  //printf("Total elements: %d \n", numElements);
  //printf("Num threads: %d \n", numThreads);
  //printf("Canonical job size: %d \n", jobSize);
  //printf("Num threads with job size augmented by one edge: %d \n", remainder);
}

void CorotationalLinearFEMMT::ComputeForceAndStiffnessMatrix(double * u, double * f, SparseMatrix * stiffnessMatrix, int warp)
{
  // launch threads
  struct CorotationalLinearFEMMT_threadArg * threadArgv = (struct CorotationalLinearFEMMT_threadArg*) malloc (sizeof(struct CorotationalLinearFEMMT_threadArg) * numThreads);

  pthread_t * tid = (pthread_t*) malloc (sizeof(pthread_t) * numThreads);

  int numVertices3 = 3 * tetMesh->getNumVertices();

  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].corotationalLinearFEMMT = this;
    threadArgv[i].u = u;
    threadArgv[i].f = &internalForceBuffer[i * numVertices3];
    threadArgv[i].stiffnessMatrix = stiffnessMatrixBuffer[i];
    threadArgv[i].warp = warp;
    threadArgv[i].rank = i;
  }

  for(int i=0; i<numThreads; i++)  
  {
    if (pthread_create(&tid[i], NULL, CorotationalLinearFEMMT_WorkerThread, &threadArgv[i]) != 0)
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

  if (f != NULL) 
    memset(f, 0, sizeof(double) * numVertices3);

  if (stiffnessMatrix != NULL) 
  {
    stiffnessMatrix->ResetToZero();

    for(int i=0; i<numThreads; i++)
    {
       double * source = &internalForceBuffer[i * numVertices3];
       if (f != NULL) 
       {
         for(int j=0; j<numVertices3; j++)
           f[j] += source[j];
       }

       *stiffnessMatrix += *(stiffnessMatrixBuffer[i]);
    }
  }
}

int CorotationalLinearFEMMT::GetStartElement(int rank)
{
  return startElement[rank];
}

int CorotationalLinearFEMMT::GetEndElement(int rank)
{
  return endElement[rank];
}

