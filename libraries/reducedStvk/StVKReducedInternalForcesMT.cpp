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

#include <pthread.h>
#include "lapack-headers.h"
#include "matrixIO.h"
#include "matrixMacros.h"
#include "matrixProjection.h"
#include "StVKReducedInternalForcesMT.h"

struct StVKReducedInternalForcesMT_threadArg
{
  StVKReducedInternalForcesMT * stVKReducedInternalForcesMT;
  double * targetBuffer;
  int rank;
};

void * StVKReducedInternalForcesMT_WorkerThread(void * arg)
{
  struct StVKReducedInternalForcesMT_threadArg * threadArgp = (struct StVKReducedInternalForcesMT_threadArg*) arg;
  StVKReducedInternalForcesMT * stVKReducedInternalForcesMT = threadArgp->stVKReducedInternalForcesMT;
  double * targetBuffer = threadArgp->targetBuffer;
  int rank = threadArgp->rank;
  int r = stVKReducedInternalForcesMT->Getr();
  int startElement = stVKReducedInternalForcesMT->GetStartElement(rank);
  int endElement = stVKReducedInternalForcesMT->GetEndElement(rank);

  double * target[3] = { &targetBuffer[0], &targetBuffer[r*StVKReducedInternalForces::GetLinearSize(r)], &targetBuffer[r*(StVKReducedInternalForces::GetLinearSize(r) + StVKReducedInternalForces::GetQuadraticSize(r))] };
  stVKReducedInternalForcesMT->ProcessElements(startElement, endElement, target);

  return NULL;
}

StVKReducedInternalForcesMT::StVKReducedInternalForcesMT(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, bool addGravity_, double g_, int numThreads_, int verbose_): 
  StVKReducedInternalForces(r, U, volumetricMesh, precomputedABCDIntegrals, 1, addGravity_, g_, verbose_), numThreads(numThreads_)
{
  internalForceBuffer = (double*) malloc (sizeof(double) * numThreads * r * (linearSize + quadraticSize + cubicSize));

  // split the workload
  int numElements = volumetricMesh->getNumElements();
  startElement = (int*) malloc (sizeof(int) * numThreads);
  endElement = (int*) malloc (sizeof(int) * numThreads);

  int remainder = numElements % numThreads;
  // the first 'remainder' nodes will process one element more
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

  if (verbose)
  {
    printf("Reduced StVK coefficients computation info:\n");
    printf("Total elements: %d \n", numElements);
    printf("Num threads: %d \n", numThreads);
    printf("Canonical job size: %d \n", jobSize);
    printf("Num threads with job size augmented by one element: %d \n", remainder);
  }

  // launch threads 
    
  struct StVKReducedInternalForcesMT_threadArg * threadArgv = (struct StVKReducedInternalForcesMT_threadArg*) malloc (sizeof(struct StVKReducedInternalForcesMT_threadArg) * numThreads);

  pthread_t * tid = (pthread_t*) malloc (sizeof(pthread_t) * numThreads);

  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].stVKReducedInternalForcesMT = this;
    threadArgv[i].targetBuffer = &internalForceBuffer[i * r * (linearSize + quadraticSize + cubicSize)];
    threadArgv[i].rank = i;
  }
    
  memset(internalForceBuffer, 0, sizeof(double) * numThreads * r * (linearSize + quadraticSize + cubicSize));

  for(int i=0; i<numThreads; i++)
  {
    if (pthread_create(&tid[i], NULL, StVKReducedInternalForcesMT_WorkerThread, &threadArgv[i]) != 0)
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

  // assemble
  memset(linearCoef_, 0, sizeof(double) * r * linearSize);
  memset(quadraticCoef_, 0, sizeof(double) * r * quadraticSize);
  memset(cubicCoef_, 0, sizeof(double) * r * cubicSize);
  for(int i=0; i<numThreads; i++)
  {
    double * sourceLinear = &internalForceBuffer[i * r * (linearSize + quadraticSize + cubicSize)];
    for(int j=0; j<r*linearSize; j++)
      linearCoef_[j] += sourceLinear[j];

    double * sourceQuadratic = &internalForceBuffer[i * r * (linearSize + quadraticSize + cubicSize) + r * linearSize];
    for(int j=0; j<r*quadraticSize; j++)
      quadraticCoef_[j] += sourceQuadratic[j];

    double * sourceCubic = &internalForceBuffer[i * r * (linearSize + quadraticSize + cubicSize) + r * (linearSize + quadraticSize)];
    for(int j=0; j<r*cubicSize; j++)
      cubicCoef_[j] += sourceCubic[j];
  }
}

StVKReducedInternalForcesMT::~StVKReducedInternalForcesMT()
{
  free(startElement);
  free(endElement);
  free(internalForceBuffer);
}

int StVKReducedInternalForcesMT::GetStartElement(int rank)
{
  return startElement[rank];
}

int StVKReducedInternalForcesMT::GetEndElement(int rank)
{
  return endElement[rank];
}


