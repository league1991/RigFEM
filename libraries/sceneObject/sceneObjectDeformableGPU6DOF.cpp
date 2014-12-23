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

#include "sceneObjectDeformableGPU6DOF.h"

SceneObjectDeformableGPU6DOF::SceneObjectDeformableGPU6DOF(const char * filenameOBJ): SceneObjectWithRestPosition(filenameOBJ), SceneObjectDeformable(filenameOBJ), SceneObjectDeformableGPU(filenameOBJ), SceneObjectDeformable6DOF(filenameOBJ)
{
}

SceneObjectDeformableGPU6DOF::~SceneObjectDeformableGPU6DOF() {}

void SceneObjectDeformableGPU6DOF::Render()
{
  double M[16] = {R[0], R[3], R[6], 0,
        R[1], R[4], R[7], 0,
        R[2], R[5], R[8], 0,
        centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObjectDeformableGPU::Render();
  glPopMatrix();
}

void SceneObjectDeformableGPU6DOF::RenderShadow(double ground[4], double light[4])
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    SetShadowingModelviewMatrix(ground, light);
    double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

    glMultMatrixd(M);
    render_coarseToFine->RenderShadow(0.2);
  glPopMatrix();
}

/*
void SceneObjectDeformableGPU6DOF::RenderVertices()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObjectReducedGPU::RenderVertices();
  glPopMatrix();
}

void SceneObjectDeformableGPU6DOF::RenderVertices_Selection()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObjectReducedGPU::RenderVertices_Selection();
  glPopMatrix();
}

void SceneObjectDeformableGPU6DOF::RenderEdges()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObjectReducedGPU::RenderEdges();
  glPopMatrix();
}

void SceneObjectDeformableGPU6DOF::RenderVertices(int numVertices, int * vertexList)
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObjectReducedGPU::RenderVertices(numVertices, vertexList);
  glPopMatrix();
}

int SceneObjectDeformableGPU6DOF::closestVertex(double queryPosX, double queryPosY, double queryPosZ, double * uBuffer)
{
  // transform the position to the local frame
  double queryPos[3] = { queryPosX, queryPosY, queryPosZ };
  double localQueryPos[3];
  TransformToLocal(queryPos, localQueryPos);

  return SceneObjectReducedGPU::closestVertex(localQueryPos[0], localQueryPos[1], localQueryPos[2], uBuffer);
}
*/
