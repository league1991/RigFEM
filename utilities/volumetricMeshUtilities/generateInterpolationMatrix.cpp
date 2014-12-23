/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "generateInterpolationMatrix" utility , Copyright (C) 2007 CMU,       *
 *                                             2009 MIT, 2014 USC        *
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
 * included with this utility in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "volumetricMesh.h"
#include "volumetricMeshLoader.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "generateInterpolationMatrix.h"
#include "objMesh.h"
#include "matrixIO.h"
#include "matrixMacros.h"
#include "getopts.h"

/*
  Builds an interpolation matrix A that can be used to efficiently transfer deformations from
  volumetric meshes to surface (obj) meshes:

  surface mesh deformation = A * volumetric mesh deformation
*/

int main(int argc, char ** argv)
{
  if (argc < 4)
  {
    printf("Generates the sparse matrix that interpolates given volumetric mesh datamatrix to vertices of the given mesh.\n");
    printf("Usage: %s <volumetric mesh file> <target obj file> <output sparse matrix> [-i interpolant file] [-z threshold]\n",argv[0]);
    printf("-i : use the specified interpolant (default: generate interpolant from scratch)\n");
    printf("-z : assign zero mode to vertices too far away from the volumetric mesh\n");
    return 0;
  }

  char * meshFile = argv[1];
  char * objMeshname = argv[2];
  char * outputFilename = argv[3];

  char zeroThresholdString[4096] = "__none";
  char interpolantFile[4096] = "__none";
  opt_t opttable[] =
  {
    { (char*)"i", OPTSTR, &interpolantFile },
    { (char*)"z", OPTSTR, &zeroThresholdString },
    { NULL, 0, NULL }
  };

  argv += 3;
  argc -= 3;
  int optup = getopts(argc,argv,opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %s.\n",argv[optup]);
    return 1;
  }

  double threshold;
  if (strcmp(zeroThresholdString,"__none") == 0)
    threshold = -1;
  else
    threshold = strtod(zeroThresholdString, NULL);

  VolumetricMesh * volumetricMesh = VolumetricMeshLoader::load(meshFile);

  if (volumetricMesh == NULL)
  {
    printf("Error: cannot load volumetric mesh %s.\n", meshFile);
    exit(1);
  }
 
  int n = volumetricMesh->getNumVertices();
  int nel = volumetricMesh->getNumElements();
  printf("Info on %s:\n",meshFile);
  printf("Num vertices: %d\n",n);
  printf("Num elements: %d\n",nel);

  ObjMesh * objMesh = new ObjMesh(objMeshname);
  int numInterpolationLocations = objMesh->getNumVertices();
  double * interpolationLocations = (double*) malloc (sizeof(double) * 3 * numInterpolationLocations);
  for(int i=0; i< numInterpolationLocations; i++)
  {
    Vec3d pos = objMesh->getPosition(i);
    interpolationLocations[3*i+0] = pos[0];
    interpolationLocations[3*i+1] = pos[1];
    interpolationLocations[3*i+2] = pos[2];
  }

  int * vertices;
  double * weights;
  if (strcmp(interpolantFile, "__none") == 0)
  {
    printf("Generating interpolation weights...\n");
    int numExternalVertices = volumetricMesh->generateInterpolationWeights(
      numInterpolationLocations, interpolationLocations, &vertices, &weights, threshold);
    printf("Encountered %d vertices not belonging to any voxel.\n", numExternalVertices);
  }
  else
  {
    printf("Loading interpolation weights from %s...\n", interpolantFile);
    int code = volumetricMesh->loadInterpolationWeights(interpolantFile, numInterpolationLocations, volumetricMesh->getNumElementVertices(), &vertices, &weights);
    if (code != 0)
    {
      printf("Error loading interpolation weights.\n");
      exit(1);
    }    
  }

  SparseMatrix * A;
  GenerateInterpolationMatrix::generate(numInterpolationLocations, volumetricMesh->getNumElementVertices(), vertices, weights, &A, volumetricMesh->getNumVertices());

  A->Save(outputFilename);

  return 0;
}

