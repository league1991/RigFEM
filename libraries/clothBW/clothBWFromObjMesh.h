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

/*
  Constructs a cloth model from an obj mesh.
*/

#ifndef _CLOTHBWFROMOBJMESH_H_
#define _CLOTHBWFROMOBJMESH_H_

#include "clothBW.h"
#include "objMesh.h"

class ClothBWFromObjMesh
{
public:
  
  // generate a cloth model from the given mesh (builds tensile, shear, and bending springs)
  // surface density and stiffnesses are the same for every triangle
  static int GenerateClothBW(ObjMesh * mesh, ClothBW ** clothBW, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffnessU, double bendStiffnessV, double damping, int bu = 1.0, int bv = 1.0, int addGravity=0);
  
  // NOTE: materialGroup 0 is hard-coded as "default" in ObjMesh.cpp. So if a .mtl file 
  // specifies 2 materials (with 2 'usemtl' calls), there will actually be 3 material groups.
  // As a result, the density/stiffness arrays that specify a value for each material group
  // all must contain a value at the beginning for the default material (even if no vertices
  // belong to this default group).
  
  // generate a cloth model from the given mesh (builds tensile, shear, and bending springs)
  // user passes array of doubles to specify surface densities and stiffness values for each material group
  static int GenerateClothBW(ObjMesh * mesh, ClothBW ** clothBW, int numMaterialGroups, double * surfaceDensity, double * tensileStiffness, double * shearStiffness, double * bendStiffnessU, double * bendStiffnessV, double * damping, int * bu = NULL, int * bv = NULL, int addGravity=0);
  
protected:
};

#endif

