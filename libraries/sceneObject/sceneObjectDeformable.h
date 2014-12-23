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
  A general deformable scene object.
  See also sceneObject.h . 
*/

#ifndef _SCENEOBJECTDEFORMABLE_H_
#define _SCENEOBJECTDEFORMABLE_H_

#include "lighting.h"
#include "sceneObjectWithRestPosition.h"

class SceneObjectDeformable : public virtual SceneObjectWithRestPosition
{
public:
  SceneObjectDeformable(const char * filenameOBJ);
  SceneObjectDeformable(ObjMesh * objMesh, bool deepCopy = true);
  virtual ~SceneObjectDeformable();

  // sets the current dynamic vertex positions to the rest position + specified deformation
  void SetVertexDeformations(double * u);
  void SetVertexDeformations(float * u);

  // adds deformations to current dynamicPosition of the vertices
  void AddVertexDeformations(double * u);

  void ResetDeformationToRest();

  inline void GetSingleVertexRestPosition(int vertex, double * x, double * y, double * z);
  inline void SetSingleVertexRestPosition(int vertex, double x, double y, double z);
  inline void GetSingleVertexPositionFromBuffer(int vertex, double * x, double * y, double * z);

  virtual void SetLighting(Lighting * lighting);

protected:

};

inline void SceneObjectDeformable::GetSingleVertexRestPosition(int vertex, double * x, double * y, double * z)
{
  *x = restPosition[3*vertex+0];
  *y = restPosition[3*vertex+1];
  *z = restPosition[3*vertex+2];
}

inline void SceneObjectDeformable::SetSingleVertexRestPosition(int vertex, double x, double y, double z)
{
  restPosition[3*vertex+0] = x;
  restPosition[3*vertex+1] = y;
  restPosition[3*vertex+2] = z;
}

inline void SceneObjectDeformable::GetSingleVertexPositionFromBuffer(int vertex, double * x, double * y, double * z)
{
  Vec3d pos = mesh->getPosition(vertex);
  *x = pos[0];
  *y = pos[1];
  *z = pos[2];
}

#endif

