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

#ifndef _OBJMESHGPUDEFORMER_UUQ_H_
#define _OBJMESHGPUDEFORMER_UUQ_H_

#include "objMeshGPUDeformer.h"

class ObjMeshGPUDeformer_uUq : public ObjMeshGPUDeformer
{
public:
  ObjMeshGPUDeformer_uUq();
  void Init(ObjMesh * mesh, ObjMeshRender * meshRender, int r, double * U, int renderingMode);
  virtual void Clone(ObjMeshGPUDeformer_uUq * ObjMeshGPUDeformer_uUq_source);
  virtual ~ObjMeshGPUDeformer_uUq() = 0;

  void Setqdv(double * q);
  void Setqfv(float * q);

  void Compute_uUq();

protected:

  char ** compilerOptions;
  int r;
  double * U;
  float * q;

  GLuint UTextureID;
  int UTextureSize;
  float invFragmentTextureSize;
  int numPixelsInARow;
  void InitUTexture();

  static bool hasStaticBeenInitialized;
  int InitializeCGShaders();

  virtual void * GetDerivedData() = 0;
  virtual void SetDerivedData(void * data) = 0;

  // fragment u=Uq shader 
  CGprofile Fragment_uUqProfile;
  CGprogram Fragment_uUqProgram; // note: this program is not static, so that can support scenes with variables number of r's
  CGparameter UtexParam;
  CGparameter qParam;
  CGparameter vertexDeformationTextureSizeParam;
  CGparameter SIZEParam;
  CGparameter texUDeltaParam;
  CGparameter SIZEdivRParam;
};

#endif

