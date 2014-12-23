/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  A multi-threaded version of the StVKReducedInternalForces class, using 
  the wxWidgets threading API. This class requires the wxWidgets header 
  files. 
*/

#include <vector>
#include "wx/wx.h"
#include "lapack-headers.h"
#include "matrixIO.h"
#include "matrixMacros.h"
#include "matrixProjection.h"
#include "StVKReducedInternalForcesWX.h"

struct StVKReducedInternalForcesWX_threadArg
{
  StVKReducedInternalForcesWX * stVKReducedInternalForcesWX;
  double * targetBuffer;
  int rank;
};

void * StVKReducedInternalForcesWX_WorkerThread(void * arg)
{
  struct StVKReducedInternalForcesWX_threadArg * threadArgp = (struct StVKReducedInternalForcesWX_threadArg*) arg;
  StVKReducedInternalForcesWX * stVKReducedInternalForcesWX = threadArgp->stVKReducedInternalForcesWX;
  double * targetBuffer = threadArgp->targetBuffer;
  int rank = threadArgp->rank;
  int r = stVKReducedInternalForcesWX->Getr();
  int startElement = stVKReducedInternalForcesWX->GetStartElement(rank);
  int endElement = stVKReducedInternalForcesWX->GetEndElement(rank);

  double * target[3] = { &targetBuffer[0], &targetBuffer[r*StVKReducedInternalForces::GetLinearSize(r)], &targetBuffer[r*(StVKReducedInternalForces::GetLinearSize(r) + StVKReducedInternalForces::GetQuadraticSize(r))] };
  stVKReducedInternalForcesWX->ProcessElements(startElement, endElement, target);

  return NULL;
}

class MyThread : public wxThread
{
public:
  MyThread(void * arg_) : wxThread(wxTHREAD_JOINABLE), arg(arg_) {}

protected:
  void * arg;
  virtual ExitCode Entry();
};

MyThread::ExitCode MyThread::Entry()
{
  StVKReducedInternalForcesWX_WorkerThread(arg);
  return (MyThread::ExitCode)0;
}

StVKReducedInternalForcesWX::StVKReducedInternalForcesWX(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, bool addGravity_, int numThreads): StVKReducedInternalForces(r, U, volumetricMesh, precomputedABCDIntegrals, 1, addGravity_)
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

  printf("Reduced StVK coefficients computation info:\n");
  printf("Total elements: %d \n",numElements);
  printf("Num threads: %d \n",numThreads);
  printf("Canonical job size: %d \n",jobSize);
  printf("Num threads with job size augmented by one element: %d \n",remainder);

  // launch threads 
  struct StVKReducedInternalForcesWX_threadArg * threadArgv = (struct StVKReducedInternalForcesWX_threadArg*) malloc (sizeof(struct StVKReducedInternalForcesWX_threadArg) * numThreads);

  for(int i=0; i<numThreads; i++)
  {
    threadArgv[i].stVKReducedInternalForcesWX = this;
    threadArgv[i].targetBuffer = &internalForceBuffer[i * r * (linearSize + quadraticSize + cubicSize)];
    threadArgv[i].rank = i;
  }
    
  memset(internalForceBuffer, 0, sizeof(double) * numThreads * r * (linearSize + quadraticSize + cubicSize));

  MyThread ** threads = (MyThread**) malloc (sizeof(MyThread*) * numThreads);

  for(int i=0; i<numThreads; i++)
  {
    threads[i] = new MyThread(&threadArgv[i]);

    if (threads[i]->Create() != wxTHREAD_NO_ERROR)
    {
      printf("Error: unable to create thread %d.\n", i);
      throw 1;
    }

    threads[i]->SetPriority(WXTHREAD_MAX_PRIORITY);

    if (threads[i]->Run() != wxTHREAD_NO_ERROR)
    {
      printf("Error: unable to run thread %d.\n", i);
      throw 2;
    }
  }

  for(int i=0; i<numThreads; i++)
  {
    if (threads[i]->Wait() != 0)
    {
      printf("Error: unable to join thread %d.\n", i);
      throw 3;
    } 

    delete(threads[i]);
  }

  free(threads);
  free(threadArgv);

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

StVKReducedInternalForcesWX::~StVKReducedInternalForcesWX()
{
  free(startElement);
  free(endElement);
  free(internalForceBuffer);
}

int StVKReducedInternalForcesWX::GetStartElement(int rank)
{
  return startElement[rank];
}

int StVKReducedInternalForcesWX::GetEndElement(int rank)
{
  return endElement[rank];
}

