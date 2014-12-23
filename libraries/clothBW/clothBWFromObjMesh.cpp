/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "clothBW" library , Copyright (C) 2014 USC                            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Andy Pierce, Yu Yu Xu, Jernej Barbic                     *
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

#include "macros.h"
#include "clothBWFromObjMesh.h"
using namespace std;

int ClothBWFromObjMesh::GenerateClothBW(ObjMesh * mesh, ClothBW ** clothBW, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffnessU, double bendStiffnessV, double damping, int bu, int bv, int addGravity)
{
  int numMaterialGroups = mesh->getNumMaterials();
  
  // create arrays of density & stiffness values
  double * surfaceDensityArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * tensileStiffnessArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * shearStiffnessArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * bendStiffnessUArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * bendStiffnessVArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  double * dampingArray = (double*) malloc (sizeof(double) * numMaterialGroups);
  int * buArray = (int*) malloc (sizeof(int) * numMaterialGroups);
  int * bvArray = (int*) malloc (sizeof(int) * numMaterialGroups);
  
  // fill arrays all with the same value
  for (int i=0; i<numMaterialGroups; i++)
  {
    surfaceDensityArray[i] = surfaceDensity;
    tensileStiffnessArray[i] = tensileStiffness;
    shearStiffnessArray[i] = shearStiffness;
    bendStiffnessUArray[i] = bendStiffnessU;
    bendStiffnessVArray[i] = bendStiffnessV;
    dampingArray[i] = damping;
    buArray[i] = bu;
    bvArray[i] = bv;
  }
  
  // construct clothBW
  GenerateClothBW(mesh, clothBW, numMaterialGroups, surfaceDensityArray, tensileStiffnessArray, shearStiffnessArray, bendStiffnessUArray, bendStiffnessVArray, dampingArray, buArray, bvArray, addGravity);
  
  // free memory we allocated
  free(surfaceDensityArray);
  free(tensileStiffnessArray);
  free(shearStiffnessArray);
  free(bendStiffnessUArray);
  free(bendStiffnessVArray);
  free(dampingArray);
  free(buArray);
  free(bvArray);

  return 0;
}

int ClothBWFromObjMesh::GenerateClothBW(ObjMesh * mesh, ClothBW ** clothBW, int numMaterialGroups, double * surfaceDensity, double * tensileStiffness, double * shearStiffness, double * bendStiffnessU, double * bendStiffnessV, double * damping, int * bu, int * bv, int addGravity)
{
  mesh->triangulate();
  
  int numParticles;
  double * restPositions;
  int numTriangles;
  int * triangles;
  int numGroups;
  int * triangleGroups;
  
  mesh->exportGeometry(&numParticles, &restPositions, &numTriangles, &triangles, &numGroups, &triangleGroups);

  if (numGroups != numMaterialGroups)
  {
    printf("Mismatch in the number of groups. Mesh has %d groups.\n", numGroups);
    return 1;
  }
  
  // compute masses
  vector<double> groupSurfaceMassDensity;
  for(int i=0; i<numGroups; i++)
    groupSurfaceMassDensity.push_back(surfaceDensity[i]);

  vector<double> massesV;
  mesh->computeMassPerVertex(groupSurfaceMassDensity, massesV);

  double * masses = (double*) malloc (sizeof(double) * numParticles);
  for(int i=0; i<numParticles; i++)
    masses[i] = massesV[i];
  
  *clothBW = new ClothBW(numParticles,  masses, restPositions, numTriangles,
                         triangles, triangleGroups, numGroups, tensileStiffness,
                         shearStiffness, bendStiffnessU, bendStiffnessV,
                         damping, addGravity);
  
  free(restPositions);
  free(triangles);
  free(triangleGroups);
  free(masses);
  
  return 0;
}

