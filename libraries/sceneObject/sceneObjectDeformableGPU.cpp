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

#include "sceneObjectDeformableGPU.h"

SceneObjectDeformableGPU::SceneObjectDeformableGPU(const char * filenameOBJ) : SceneObjectWithRestPosition(filenameOBJ), SceneObjectDeformable(filenameOBJ)
{
  try
  {
    render_coarseToFine = new ObjMeshGPUDeformer_coarseToFine_fbo();
    //SetMaterialAlpha(0.5);
  }
  catch(int exception)
  {
    printf("Failed to initialize the interp GP-GPU renderer. Exception code: %d\n", exception);
    throw 1;
  }
}

void SceneObjectDeformableGPU::SetInterpolation(int numCoarseVertices, int numElementVertices, int * vertices, double * weights)
{
  try
  {
    render_coarseToFine->Init(mesh, meshRender, numCoarseVertices, numElementVertices, vertices, weights, renderMode);
  }
  catch(int exception)
  {
    printf("Failed to initialize the interp GP-GPU renderer (interp). Exception code: %d\n", exception);
    throw 1;
  }
}

SceneObjectDeformableGPU::~SceneObjectDeformableGPU()
{
  delete(render_coarseToFine);
}

void SceneObjectDeformableGPU::SetCoarseDeformations(double * u)
{
  render_coarseToFine->SetCoarseDeformations(u);
}

void SceneObjectDeformableGPU::Render()
{
  render_coarseToFine->Render();
}

void SceneObjectDeformableGPU::RenderShadow(double ground[4], double light[4])
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  SetShadowingModelviewMatrix(ground, light);
  render_coarseToFine->RenderShadow(0.2);
  glPopMatrix();
}

void SceneObjectDeformableGPU::SetLighting(Lighting * lighting)
{
  //setGPULighting.SetLighting(render_coarseToFine, lighting);
}

void SceneObjectDeformableGPU::Getu(double * u)
{
  render_coarseToFine->ReadBack_u(u);
}

