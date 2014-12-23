/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "objMeshGPUDeformer" library , Copyright (C) 2007 CMU, 2009 MIT,      *
 *                                                        2014 USC       *
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

#ifndef _OBJMESHGPUDEFORMER_COARSETOFINE_H_
#define _OBJMESHGPUDEFORMER_COARSETOFINE_H_

#include "objMeshGPUDeformer.h"

class ObjMeshGPUDeformer_coarseToFine : public ObjMeshGPUDeformer
{
public:
  ObjMeshGPUDeformer_coarseToFine();
  void Init(ObjMesh * mesh, ObjMeshRender * meshRender, int numCoarseVertices, int interp_numElementVertices, int * interp_vertices, double * interp_weights, int renderingMode);
  virtual ~ObjMeshGPUDeformer_coarseToFine() = 0;

  void SetCoarseDeformations(double * uCoarse);
  void PrintFineDeformationTexture();
  void PrintCoarseDeformationTexture();

protected:

  int interp_numElementVertices; 
  int * interp_vertices; 
  double * interp_weights;
  float * uCoarseTextureData;

  int numCoarseVertices;

  GLuint interpolationTextureID;
  GLuint coarseDeformationTextureID;

  int coarseDeformationTextureSize;
  int interpolationTextureSize;
  float invInterpolationTextureSize;

  void InitInterpolationTexture();
  void InitCoarseDeformationTexture();

  void InterpolateDeformations();
  void UploadDeformations();

  static bool hasStaticBeenInitialized;
  int InitializeCGShaders();

  virtual void * GetDerivedData() = 0;
  virtual void SetDerivedData(void * data) = 0;

  // fragment shaders
  static CGprofile Fragment_InterpolationProfile;
  static CGprogram Fragment_InterpolationProgram; 
  static CGprogram FragmentPass2ProgramNoTexture;

  // fragment shader parameters
  CGparameter interpolationTextureParam;
  CGparameter interpolationTextureSizeParam;
  CGparameter coarseDeformationTextureParam;
  CGparameter coarseDeformationTextureSizeParam;
  CGparameter vertexDeformationTextureSizeParam;

  int internalFormat;
};

#endif

