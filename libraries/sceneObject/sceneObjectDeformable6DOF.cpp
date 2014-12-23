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
#include "sceneObjectDeformable6DOF.h"

SceneObjectDeformable6DOF::SceneObjectDeformable6DOF(const char * filenameOBJ): SceneObjectWithRestPosition(filenameOBJ), SceneObjectDeformable(filenameOBJ), SceneObject6DOF(filenameOBJ)
{
}

SceneObjectDeformable6DOF::SceneObjectDeformable6DOF(ObjMesh * objMesh, bool deepCopy): SceneObjectWithRestPosition(objMesh, deepCopy), SceneObjectDeformable(objMesh, deepCopy), SceneObject6DOF(objMesh, deepCopy)
{
}

void SceneObjectDeformable6DOF::GetSingleVertexPosition(int vertex, double * x, double * y, double * z)
{
  Vec3d pos = mesh->getPosition(vertex);
  double x0 = pos[0];
  double y0 = pos[1];
  double z0 = pos[2];

  // transform x0, y0, z0 via centerOfMass, R
  // x = centerOfMass + R * x0
  *x = R[0] * x0 + R[1] * y0 + R[2] * z0 + centerOfMass[0];
  *y = R[3] * x0 + R[4] * y0 + R[5] * z0 + centerOfMass[1];
  *z = R[6] * x0 + R[7] * y0 + R[8] * z0 + centerOfMass[2];
}

SceneObjectDeformable6DOF::~SceneObjectDeformable6DOF() {}

void SceneObjectDeformable6DOF::Render() 
{ 
  SceneObject6DOF::Render(); 
}

int SceneObjectDeformable6DOF::GetClosestVertex(Vec3d & queryPos, double * distance, double * auxVertexBuffer)
{
  // transform the position to the local frame
  double queryPosv[3] = { queryPos[0], queryPos[1], queryPos[2] };
  double localQueryPosv[3];
  TransformToLocal(queryPosv, localQueryPosv);
  Vec3d localQueryPos(localQueryPosv);

  return SceneObjectDeformable::GetClosestVertex(localQueryPos, distance, auxVertexBuffer);
}

