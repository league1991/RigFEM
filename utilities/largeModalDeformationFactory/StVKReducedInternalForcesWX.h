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

#ifndef _STVKREDUCEDINTERNALFORCESWX_H_
#define _STVKREDUCEDINTERNALFORCESWX_H_

#include "StVKReducedInternalForces.h"

class StVKReducedInternalForcesWX : public StVKReducedInternalForces
{
public:
  // creates the reduced coefficients
  StVKReducedInternalForcesWX(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, bool addGravity, int numThreads);

  ~StVKReducedInternalForcesWX();

  int GetStartElement(int rank);  
  int GetEndElement(int rank);

protected:
  int numThreads;
  int * startElement, * endElement;
  double * internalForceBuffer;
};

#endif

