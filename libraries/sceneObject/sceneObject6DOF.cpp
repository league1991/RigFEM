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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "sceneObject6DOF.h"

SceneObject6DOF::SceneObject6DOF(const char * filenameOBJ): SceneObjectWithRestPosition(filenameOBJ), scale(1.0) 
{
  Construct();
}

SceneObject6DOF::SceneObject6DOF(ObjMesh * objMesh, bool deepCopy): SceneObjectWithRestPosition(objMesh, deepCopy), scale(1.0) 
{
  Construct();
}

void SceneObject6DOF::Construct()
{
  memset(centerOfMass, 0, sizeof(double) * 3);
  memset(R, 0, sizeof(double) * 9);
  R[0] = 1; R[4] = 1; R[8] = 1;
  rotation = 1;
}

SceneObject6DOF::~SceneObject6DOF()
{
}

void SceneObject6DOF::SetScale(double scale)
{
  this->scale = scale;
}

void SceneObject6DOF::GetSingleVertexPosition(int vertex, double * x, double * y, double * z)
{
  double x0 = restPosition[3*vertex+0];
  double y0 = restPosition[3*vertex+1];
  double z0 = restPosition[3*vertex+2];

  // transform x0, y0, z0 via centerOfMass, R
  // x = centerOfMass + R * x0
  *x = R[0] * x0 + R[1] * y0 + R[2] * z0 + centerOfMass[0];
  *y = R[3] * x0 + R[4] * y0 + R[5] * z0 + centerOfMass[1];
  *z = R[6] * x0 + R[7] * y0 + R[8] * z0 + centerOfMass[2];
}

void SceneObject6DOF::GetSingleVertexVelocity(int vertex, double * objectVel, double * objectAngVel, double * velx, double * vely, double * velz)
{
  *velx = objectVel[0];
  *vely = objectVel[1];
  *velz = objectVel[2];

  // the handle (vector from the center of mass to the point, in the world-coordinate frame)
  double x0 = restPosition[3*vertex+0];
  double y0 = restPosition[3*vertex+1];
  double z0 = restPosition[3*vertex+2];
  double hx0 = R[0] * x0 + R[1] * y0 + R[2] * z0; 
  double hy0 = R[3] * x0 + R[4] * y0 + R[5] * z0;
  double hz0 = R[6] * x0 + R[7] * y0 + R[8] * z0;

  CROSSPRODUCT_ADD(objectAngVel[0], objectAngVel[1], objectAngVel[2], hx0, hy0, hz0, *velx, *vely, *velz);
}

void SceneObject6DOF::SetOpenGLModelviewMatrix()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
          R[2], R[5], R[8], 0,
          centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMultMatrixd(M);
}

void SceneObject6DOF::Render()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glMultMatrixd(M);
    SceneObject::Render();
  glPopMatrix();
}

void SceneObject6DOF::RenderShadow(double ground[4], double light[4])
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  SetShadowingModelviewMatrix(ground, light);

  double M[16] = {R[0], R[3], R[6], 0,
                  R[1], R[4], R[7], 0,
                  R[2], R[5], R[8], 0,
                  centerOfMass[0], centerOfMass[1], centerOfMass[2], 1 };

  glMultMatrixd(M);
  
  bool texEnabled = AreTexturesEnabled();
  DisableTextures();
  SceneObject::Render();
  if (texEnabled)
    EnableTextures();

  glPopMatrix();
}

void SceneObject6DOF::RenderVertices()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObject::RenderVertices();
  glPopMatrix();
}

void SceneObject6DOF::RenderVertices_Selection()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObject::RenderVertices_Selection();
  glPopMatrix();
}

void SceneObject6DOF::RenderEdges()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObject::RenderEdges();
  glPopMatrix();
}

void SceneObject6DOF::RenderNormals()
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrixd(M);
    SceneObject::RenderNormals();
  glPopMatrix();
}

void SceneObject6DOF::RenderVertices(int numVertices, int * vertexList)
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glMultMatrixd(M);
	SceneObject::RenderVertices(numVertices, vertexList);
  glPopMatrix();
}

void SceneObject6DOF::RenderLocalFrame(double axesSize)
{
  double M[16] = {R[0], R[3], R[6], 0,
          R[1], R[4], R[7], 0,
	  R[2], R[5], R[8], 0,
	  centerOfMass[0],    centerOfMass[1],    centerOfMass[2], 1 };

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glMultMatrixd(M);

    glScalef(axesSize,axesSize,axesSize);
    glBegin(GL_LINES);
    glColor3f(1,0,0);
      glVertex3f(0,0,0);
      glVertex3f(1,0,0);
    glColor3f(0,1,0);
      glVertex3f(0,0,0);
      glVertex3f(0,1,0);
    glColor3f(0,0,1);
      glVertex3f(0,0,0);
      glVertex3f(0,0,1);
    glEnd();

  glPopMatrix();
}

void SceneObject6DOF::TransformToLocal(double * globalVector, double * localVector)
{
  double temp[3] = { globalVector[0] - centerOfMass[0], globalVector[1] - centerOfMass[1], globalVector[2] - centerOfMass[2] };

  if (rotation != 0)
  {
    // global = R * local + pos
    // local = R^{-1} * (global - pos)

    // we know that R = s * Rot, so we can invert R as follows:  
    // R^{-1} = 1 / s Rot^T = 1 / s^2 * R

    localVector[0] = R[0] * temp[0] + R[3] * temp[1] + R[6] * temp[2];
    localVector[1] = R[1] * temp[0] + R[4] * temp[1] + R[7] * temp[2];
    localVector[2] = R[2] * temp[0] + R[5] * temp[1] + R[8] * temp[2];
    
    double invScale2 = 1.0 / (scale * scale);
    localVector[0] *= invScale2;
    localVector[1] *= invScale2;
    localVector[2] *= invScale2;
  }
  else
  {
    // global = R * local + pos
    // local = R^{-1} * (global - pos)
    // R is not necessarily a rotation

    double RInv[9];
    inverse3x3(R, RInv);

    localVector[0] = RInv[0] * temp[0] + RInv[1] * temp[1] + RInv[2] * temp[2];
    localVector[1] = RInv[3] * temp[0] + RInv[4] * temp[1] + RInv[5] * temp[2];
    localVector[2] = RInv[6] * temp[0] + RInv[7] * temp[1] + RInv[8] * temp[2];
  }
}

void SceneObject6DOF::TransformToGlobal(double * localVector, double * globalVector)
{
  globalVector[0] = centerOfMass[0] + R[0] * localVector[0] + R[1] * localVector[1] + R[2] * localVector[2];
  globalVector[1] = centerOfMass[1] + R[3] * localVector[0] + R[4] * localVector[1] + R[5] * localVector[2];
  globalVector[2] = centerOfMass[2] + R[6] * localVector[0] + R[7] * localVector[1] + R[8] * localVector[2];
}

int SceneObject6DOF::GetClosestVertex(Vec3d & queryPos, double * distance, double * auxVertexBuffer)
{
  // transform the position to the local frame
  double queryPosv[3] = { queryPos[0], queryPos[1], queryPos[2] };
  double localQueryPosv[3];
  TransformToLocal(queryPosv, localQueryPosv);
  Vec3d localQueryPos(localQueryPosv);

  return SceneObject::GetClosestVertex(localQueryPos, distance, auxVertexBuffer);
}

