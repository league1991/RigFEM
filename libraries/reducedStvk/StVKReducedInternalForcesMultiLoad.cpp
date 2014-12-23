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

#include "StVKReducedInternalForcesMultiLoad.h"

int StVKReducedInternalForcesMultiLoad(const char * filename, int * numModels, StVKReducedInternalForces *** stVKReducedInternalForces, int verbose)
{
  FILE * fin = fopen(filename, "rb");
  if (!fin)
  {
    printf("Error: unable to access file %s.\n", filename);
    return 1;
  }

  if ((int)(fread(numModels,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from the input cubic polynomial multifile.\n");
    return 1;
  }

  *stVKReducedInternalForces = (StVKReducedInternalForces**) malloc (sizeof(StVKReducedInternalForces*) * *numModels);
  for(int i=0; i<*numModels; i++)
  {
    int rTarget = -1;
    int bigEndianMachine = 0;
    (*stVKReducedInternalForces)[i] = new StVKReducedInternalForces(fin, rTarget, bigEndianMachine, verbose);
  }

  fclose(fin);

  return 0;
}

int StVKReducedInternalForcesMultiSave(const char * filename, int numModels, StVKReducedInternalForces ** stVKReducedInternalForces, int verbose)
{
  FILE * fout = fopen(filename, "wb");
  if (!fout)
  {
    printf("Error: unable to access file %s.\n", filename);
    return 1;
  }

  if ((int)(fwrite(&numModels,sizeof(int),1,fout)) < 1)
    return 1;
  
  for(int i=0; i<numModels; i++)
  {
    if (stVKReducedInternalForces[i] != NULL)
      stVKReducedInternalForces[i]->Save(fout);
    else
      StVKReducedInternalForces::SaveEmptyCub(fout);
  }

  fclose(fout);

  return 0;
}

