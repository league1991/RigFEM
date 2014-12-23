/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "sceneObject" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
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
  A general deformable scene object, rendered on GPU via interpolation
  See also sceneObject.h . 
*/

#ifndef _SCENEOBJECTDEFORMABLEGPU_H_
#define _SCENEOBJECTDEFORMABLEGPU_H_

#define GL_GLEXT_PROTOTYPES 1

#include "objMeshGPUDeformer_coarseToFine_fbo.h"
#include "lighting.h"
#include "sceneObjectDeformable.h"

class SceneObjectDeformableGPU : public virtual SceneObjectDeformable
{
public:
  SceneObjectDeformableGPU(const char * filenameOBJ);
  virtual ~SceneObjectDeformableGPU();

  void SetInterpolation(int numCoarseVertices, int numElementVertices, int * vertices, double * weights);

  void SetCoarseDeformations(double * u); 

  virtual void Render();
  virtual void RenderShadow(double ground[4], double light[4]);
  //virtual void RenderVertices();
  //virtual void RenderEdges();

  virtual void SetLighting(Lighting * lighting);

  virtual void Getu(double * u);

protected:
  ObjMeshGPUDeformer_coarseToFine_fbo * render_coarseToFine;
};

#endif

