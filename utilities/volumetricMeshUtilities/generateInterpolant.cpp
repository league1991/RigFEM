/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "generateInterpolant" utility , Copyright (C) 2007 CMU, 2009 MIT,     *
 *                                               2014 USC                *
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
#include <vector>
using namespace std;

#include "float.h"
#include "volumetricMesh.h"
#include "volumetricMeshLoader.h"
#include "tetMesh.h"
#include "cubicMesh.h"
#include "objMesh.h"
#include "matrixIO.h"
#include "matrixMacros.h"
#include "getopts.h"
#include "loadList.h"

/*
  Builds an interpolant that can be used to efficiently transfer deformations from
  volumetric meshes to surface (obj) meshes.
*/

int main(int argc, char ** argv)
{
  if (argc < 4)
  {
    printf("Generates an interpolant between a given volumetric mesh and a surface obj mesh.\n");
    printf("Usage: %s <volumetric mesh file> <target obj file> <output interpolant file> [-s volumetric mesh element list output filename] [-z threshold] [-S] [-T]\n",argv[0]);
    printf("-s : output list of (1-indexed) volumetric mesh elements that contain at least one obj mesh vertex\n");
    printf("-z : assign zero interpolation to vertices too far away from the volumetric mesh\n");
    return 0;
  }

  char * meshFile = argv[1];
  char * objMeshname = argv[2];
  char * outputFilename = argv[3];

  char outputElementFilename[4096] = "__none";
  char zeroThresholdString[4096] = "__none";

  opt_t opttable[] =
  {
    { (char*)"s", OPTSTR, &outputElementFilename },
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
    printf("Error: unable to load the volumetric mesh from %s.\n", meshFile);
    exit(1);
  }
 
  int n = volumetricMesh->getNumVertices();
  int nel = volumetricMesh->getNumElements();
  printf("Info on %s:\n", meshFile);
  printf("Num vertices: %d\n", n);
  printf("Num elements: %d\n", nel);

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

  int * elementList;
  int ** elementListp = NULL;
  if (strcmp(outputElementFilename, "__none") != 0)
  {
    elementListp = &elementList;    
  }

  int verbose = 1;
  int numExternalVertices;
  numExternalVertices = volumetricMesh->generateInterpolationWeights(numInterpolationLocations, interpolationLocations, &vertices, &weights, threshold, elementListp, verbose);

  printf("Saving weights to %s...\n", outputFilename); fflush(NULL);
  volumetricMesh->saveInterpolationWeights(outputFilename, numInterpolationLocations, volumetricMesh->getNumElementVertices(), vertices, weights);

  if (strcmp(outputElementFilename, "__none") != 0)
  {
    set<int> uniqueElementSet;
    for(unsigned int i=0; i<numInterpolationLocations; i++)
      uniqueElementSet.insert(elementList[i]);

    vector<int> uniqueElementList;
    for(set<int>::iterator iter = uniqueElementSet.begin(); iter != uniqueElementSet.end(); iter++)
      uniqueElementList.push_back(*iter);

    LoadList saveList;
    sort(uniqueElementList.begin(), uniqueElementList.end());
    int oneIndexed = 1;
    saveList.save(outputElementFilename, uniqueElementList.size(), &uniqueElementList[0], oneIndexed);
  }
  printf("\n");

  return 0;
}

