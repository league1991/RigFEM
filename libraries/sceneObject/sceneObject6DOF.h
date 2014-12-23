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
  A scene object that can undergo a 6-DOF rigid 3D transformation.
  See also sceneObject.h . 
*/

#ifndef _SCENEOBJECT6DOF_H_
#define _SCENEOBJECT6DOF_H_

#include "sceneObjectWithRestPosition.h"

class SceneObject6DOF : public virtual SceneObjectWithRestPosition
{
public:

  SceneObject6DOF(const char * filenameOBJ); 
  SceneObject6DOF(ObjMesh * objMesh, bool deepCopy = true); 
  virtual ~SceneObject6DOF();

  inline void SetRigidBodyParameters(double * centerOfMass, double * R);
  inline void SetRigidBodyParameters(float * centerOfMass, float * R);
  inline void GetRigidBodyParameters(double * centerOfMass, double * R);

  virtual void GetSingleVertexPosition(int vertex, double * x, double * y, double * z);
  virtual void GetSingleVertexVelocity(int vertex, double * objectVel, double * objectAngVel, double * velx, double * vely, double * velz);

  virtual void Render();
  virtual void RenderVertices();
  virtual void RenderVertices_Selection();
  virtual void RenderEdges();
  virtual void RenderVertices(int numVertices, int * vertexList);
  virtual void RenderShadow(double ground[4], double light[4]);
  void RenderLocalFrame(double axesSize = 1.0);
  void RenderNormals();

  void SetOpenGLModelviewMatrix(); // sets the OpenGL modelview matrix to the current frame of the object

  void TransformToLocal(double * globalVector, double * localVector); // transform a world-coordinate vector to the frame of the rigid object
  void TransformToGlobal(double * localVector, double * globalVector); // transform a vector from the frame of the rigid object to global world-coordinate
  inline void SetRotationMode(int rotation); // to support non-rigid transformations in "TransformToLocal"; default: rotatation=1; call this function with rotation=0 if the matrix R is not a rotation, in order to get the correct results for "TransformToLocal"
  void SetScale(double scale); // to accelerate "TransformToLocal" for transformations of the form R = scale * Rot, where Rot is a rotation; if R is in this format, you can set the scale (and set rotation=1 via "SetRotationMode" above), which increases the speed of "TransformToLocal"
  double GetScale() const { return scale; }

  int GetClosestVertex(Vec3d & queryPos, double * distance, double * auxVertexBuffer);

protected:
  void Construct();
  double centerOfMass[3];
  double R[9];
  double scale;
  int rotation;
};

inline void SceneObject6DOF::SetRigidBodyParameters(double * centerOfMass_, double * R_)
{
  centerOfMass[0] = centerOfMass_[0]; 
  centerOfMass[1] = centerOfMass_[1]; 
  centerOfMass[2] = centerOfMass_[2]; 

  R[0] = R_[0];
  R[1] = R_[1];
  R[2] = R_[2];
  R[3] = R_[3];
  R[4] = R_[4];
  R[5] = R_[5];
  R[6] = R_[6];
  R[7] = R_[7];
  R[8] = R_[8];
}

inline void SceneObject6DOF::SetRigidBodyParameters(float * centerOfMass_, float * R_)
{
  centerOfMass[0] = centerOfMass_[0]; 
  centerOfMass[1] = centerOfMass_[1]; 
  centerOfMass[2] = centerOfMass_[2]; 

  R[0] = R_[0];
  R[1] = R_[1];
  R[2] = R_[2];
  R[3] = R_[3];
  R[4] = R_[4];
  R[5] = R_[5];
  R[6] = R_[6];
  R[7] = R_[7];
  R[8] = R_[8];
}


inline void SceneObject6DOF::GetRigidBodyParameters(double * centerOfMass_, double * R_)
{
  centerOfMass_[0] = centerOfMass[0]; 
  centerOfMass_[1] = centerOfMass[1]; 
  centerOfMass_[2] = centerOfMass[2]; 

  R_[0] = R[0];
  R_[1] = R[1];
  R_[2] = R[2];
  R_[3] = R[3];
  R_[4] = R[4];
  R_[5] = R[5];
  R_[6] = R[6];
  R_[7] = R[7];
  R_[8] = R[8];
}

inline void SceneObject6DOF::SetRotationMode(int rotation_)
{
  rotation = rotation_;
}

#endif

