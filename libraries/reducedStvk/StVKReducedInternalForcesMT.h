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
  This is a multi-threaded version of the StVKReducedInternalForces class.
  It uses POSIX threads ("pthreads") as the threading API.

  The computation of the cubic polynomials is multi-threaded. Speedup is in practice very good, almost linear. We have seen computation time reductions from 30 minutes to 4 minutes when using 8 CPU cores.

  The evaluation of the reduced forces is not multi-threaded in this implementation
  (although your BLAS implementation might invoke multi-threading automatically with BLAS 2 and BLAS 3 routines).

  See also StVKReducedInternalForces.h and StVKInternalForces.h .
*/

#ifndef _STVKREDUCEDINTERNALFORCESMT_H_
#define _STVKREDUCEDINTERNALFORCESMT_H_

#include "StVKReducedInternalForces.h"

class StVKReducedInternalForcesMT : public StVKReducedInternalForces
{
public:
  // creates the reduced coefficients
  StVKReducedInternalForcesMT(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, bool addGravity, double g, int numThreads, int verbose);

  ~StVKReducedInternalForcesMT();

  int GetStartElement(int rank);  
  int GetEndElement(int rank);

protected:
  int numThreads;
  int * startElement, * endElement;
  double * internalForceBuffer;
};

#endif

