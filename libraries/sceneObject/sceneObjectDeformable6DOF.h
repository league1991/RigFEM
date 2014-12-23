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

#ifndef _SCENEOBJECTDEFORMABLE6DOF_H_
#define _SCENEOBJECTDEFORMABLE6DOF_H_

#include "sceneObject6DOF.h"
#include "sceneObjectDeformable.h"

class SceneObjectDeformable6DOF : public virtual SceneObjectDeformable, public SceneObject6DOF
{
public:

  SceneObjectDeformable6DOF(const char * filenameOBJ); 
  SceneObjectDeformable6DOF(ObjMesh * objMesh, bool deepCopy = true); 
  virtual ~SceneObjectDeformable6DOF();

  virtual void GetSingleVertexPosition(int vertex, double * x, double * y, double * z);

  virtual void Render();
  virtual void RenderVertices() { SceneObject6DOF::RenderVertices(); }
  virtual void RenderVertices_Selection() { SceneObject6DOF::RenderVertices_Selection(); }
  virtual void RenderEdges() { SceneObject6DOF::RenderEdges(); }
  virtual void RenderVertices(int numVertices, int * vertexList) { SceneObject6DOF::RenderVertices(numVertices, vertexList); }
  virtual void RenderShadow(double ground[4], double light[4]) { SceneObject6DOF::RenderShadow(ground, light); }

  int GetClosestVertex(Vec3d & queryPos, double * distance, double * auxVertexBuffer);

protected:
};

#endif

