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
  Load/save many internal force models from one file.
*/

#ifndef _STVKREDUCEDINTERNALFORCESMULTILOAD_H_
#define _STVKREDUCEDINTERNALFORCESMULTILOAD_H_

#include "StVKReducedInternalForces.h"

int StVKReducedInternalForcesMultiLoad(const char * filename, int * numModels, StVKReducedInternalForces *** stVKReducedInternalForces, int verbose=0);
int StVKReducedInternalForcesMultiSave(const char * filename, int numModels, StVKReducedInternalForces ** stVKReducedInternalForces, int verbose=0);

#endif

