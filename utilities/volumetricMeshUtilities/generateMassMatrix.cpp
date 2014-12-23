/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "generateMassMatrix" utility , Copyright (C) 2007 CMU, 2009 MIT,      *
 *                                              2014 USC                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic                                           *
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
 * included with this utility in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  Generates the mass matrix for a given volumetric mesh (cubic or tet mesh).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "generateMassMatrix.h"
#include "volumetricMeshLoader.h"

int main(int argc, char ** argv)
{
  if ( argc != 3 )
  {
    printf("Computes the mass matrix for a given volumetric mesh.\n");
    printf("Usage: %s <volumetric mesh file> <output sparse mass matrix>\n",argv[0]);
    return 1;
  }

  VolumetricMesh * volumetricMesh = VolumetricMeshLoader::load(argv[1]);
  printf("The number of vertices is: %d\n", volumetricMesh->getNumVertices());
  printf("The number of elements is: %d\n", volumetricMesh->getNumElements());

  SparseMatrix * massMatrix;
  GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, false);

  massMatrix->Save(argv[2]);

  return 0;
}

