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
#include "sceneObject.h"
#include "objMeshRender.h"
#include "objMeshEncode.h"

SceneObject::SceneObject(const char * filename):
  mesh(NULL), meshRender(NULL), displayList(0), displayListExists(false), displayListEdges(0), displayListEdgesExists(false)
{ 
  int verbose = 0;
  mesh = new ObjMesh(std::string(filename), ObjMesh::ASCII, verbose);

  int encStart = strlen(filename) - 4;
  if ((encStart > 0) && (strcmp(&filename[encStart], ".enc") == 0))
  {
    // must decode
    printf("Decoding mesh.\n");
    objMeshDecode(mesh);
    printf("Decoded mesh.\n");
  }

  Construct();
}

SceneObject::SceneObject(ObjMesh * objMesh, bool deepCopy_):
  deepCopy(deepCopy_), mesh(NULL), meshRender(NULL), displayList(0), displayListExists(false), displayListEdges(0), displayListEdgesExists(false)
{ 
  if (deepCopy)
    mesh = new ObjMesh(*objMesh);
  else
    mesh = objMesh;

  Construct();
}

void SceneObject::Construct()
{
  meshRender = new ObjMeshRender(mesh);
  if (meshRender->numTextures() > 0)
    hasTextures_ = true;
  else
    hasTextures_ = false;

  BuildFaceNormals();

  n = mesh->getNumVertices();
  renderMode = OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL;
}

SceneObject::~SceneObject()
{
  PurgeDisplayList();
  if (deepCopy)
    delete(mesh);
  delete(meshRender);
}

void SceneObject::SetMaterialAlpha(double alpha)
{
  mesh->setMaterialAlpha(alpha);
}

void SceneObject::PurgeDisplayList()
{
  if (displayListExists)
    glDeleteLists(displayList, 1);
  displayListExists = false;

  if (displayListEdgesExists)
    glDeleteLists(displayListEdges, 1);
  displayListEdgesExists = false;
}

void SceneObject::BuildDisplayList()
{
  GLenum errorCode;
  const GLubyte * errorString;
  
  errorCode = glGetError();
  if (errorCode != GL_NO_ERROR)
  {
    errorString = gluErrorString(errorCode);
    printf("OpenGL Error (start of BuildDisplayList): %s\n", errorString);
  }

  if (displayListExists)
    glDeleteLists(displayList, 1);
  
  displayList = meshRender->createDisplayList(OBJMESHRENDER_TRIANGLES, renderMode);
  displayListExists = true;

  if (displayListEdgesExists)
    glDeleteLists(displayListEdges, 1);

  displayListEdges = meshRender->createDisplayList(OBJMESHRENDER_EDGES, renderMode);
  displayListEdgesExists = true;

  errorCode = glGetError();
  if (errorCode != GL_NO_ERROR)
  {
    errorString = gluErrorString(errorCode);
    printf("OpenGL Error (end of BuildDisplayList): %s\n", errorString);
  }
}

// assumes pre-existing face normals
// second parameter is treshold angle for hard edges
void SceneObject::BuildVertexNormals(double thresholdAngle)
{
  //do stuff with structure
  mesh->buildVertexNormals(thresholdAngle);
}

void SceneObject::BuildFaceNormals()
{
  mesh->buildFaceNormals();
}

void SceneObject::BuildNormals(double thresholdAngle)
{
  BuildFaceNormals();
  BuildVertexNormals(thresholdAngle);
}

void SceneObject::SetNormalsToFaceNormals()
{
  mesh->setNormalsToFaceNormals();
}

void SceneObject::BuildNormalsFancy(double thresholdAngle)
{
  BuildFaceNormals();
  mesh->buildVertexNormalsFancy(thresholdAngle);
}

void SceneObject::Render()
{
  if(displayListExists)
    glCallList(displayList);
  else
    meshRender->render(OBJMESHRENDER_TRIANGLES, renderMode);
}

void SceneObject::SetShadowingModelviewMatrix(double ground[4], double light[4])
{
  double dot;
  double shadowMat[4][4];

  dot = ground[0] * light[0] + ground[1] * light[1] + ground[2] * light[2] + ground[3] * light[3];

  shadowMat[0][0] = dot - light[0] * ground[0];
  shadowMat[1][0] = 0.0 - light[0] * ground[1];
  shadowMat[2][0] = 0.0 - light[0] * ground[2];
  shadowMat[3][0] = 0.0 - light[0] * ground[3];

  shadowMat[0][1] = 0.0 - light[1] * ground[0];
  shadowMat[1][1] = dot - light[1] * ground[1];
  shadowMat[2][1] = 0.0 - light[1] * ground[2];
  shadowMat[3][1] = 0.0 - light[1] * ground[3];

  shadowMat[0][2] = 0.0 - light[2] * ground[0];
  shadowMat[1][2] = 0.0 - light[2] * ground[1];
  shadowMat[2][2] = dot - light[2] * ground[2];
  shadowMat[3][2] = 0.0 - light[2] * ground[3];

  shadowMat[0][3] = 0.0 - light[3] * ground[0];
  shadowMat[1][3] = 0.0 - light[3] * ground[1];
  shadowMat[2][3] = 0.0 - light[3] * ground[2];
  shadowMat[3][3] = dot - light[3] * ground[3];

  glMultMatrixd((const GLdouble*)shadowMat);
}

void SceneObject::RenderShadow(double ground[4], double light[4]) 
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  
  SetShadowingModelviewMatrix(ground, light);

  bool texEnabled = AreTexturesEnabled();
  DisableTextures();

  if(displayListExists)
    glCallList(displayList);
  else
    meshRender->render(OBJMESHRENDER_TRIANGLES, renderMode);

  if (texEnabled)
    EnableTextures();

  glPopMatrix();
}

void SceneObject::RenderVertices()
{
  meshRender->render(OBJMESHRENDER_VERTICES, renderMode);
}

void SceneObject::RenderVertices_Selection()
{
  meshRender->render(OBJMESHRENDER_VERTICES, renderMode);
}

void SceneObject::RenderEdges()
{
  if(displayListEdgesExists)
    glCallList(displayListEdges);
  else
    meshRender->render(OBJMESHRENDER_EDGES, renderMode);
}

void SceneObject::RenderFacesAndEdges()
{
  meshRender->render(OBJMESHRENDER_TRIANGLES | OBJMESHRENDER_EDGES, renderMode);
}

void SceneObject::RenderEdgesInGroup(const char * groupName)
{
  meshRender->renderGroupEdges(groupName);
}

void SceneObject::RenderVertices(int numVertices, int * vertexList)
{
  meshRender->renderSpecifiedVertices(vertexList, numVertices);
}

void SceneObject::RenderVertex(int vertex)
{
  meshRender->renderVertex(vertex);
}

// shows all point labels
void SceneObject::ShowPointLabels()
{
  ShowPointLabels(0,(int)(mesh->getNumVertices())-1);
}

// shows point labels from [k to l]
void SceneObject::ShowPointLabels(int k, int l)
{
  glColor3f(0,0,0);

  // show point labels
  // labels are printed out in the range 1... , not 0...
  for (int i=k; i<= l; i++)
  {
    Vec3d pos = mesh->getPosition(i);
    PrintBitmapInteger(pos[0], pos[1], pos[2], i+1);
  }
}

int SceneObject::GetClosestVertex(Vec3d & queryPos, double * distance, double * auxVertexBuffer)
{
  return mesh->getClosestVertex(queryPos, distance);
}

// highlights vertex i, i=0,1,2,...,n-1
void SceneObject::HighlightVertex(int i)
{
  glColor3f(0,1,0);
  glPointSize(8.0);

  Vec3d pos = mesh->getPosition(i);

  glBegin(GL_POINTS);
    glVertex3f(pos[0], pos[1], pos[2]);
  glEnd();
}

bool SceneObject::AreTexturesEnabled()
{
  return ((renderMode & OBJMESHRENDER_TEXTURE) != 0);
}

void SceneObject::EnableTextures()
{
  renderMode = renderMode | OBJMESHRENDER_TEXTURE;
}

void SceneObject::DisableTextures()
{
  renderMode = renderMode & (~OBJMESHRENDER_TEXTURE);
}

int SceneObject::SetUpTextures(LightingModulationType lightingModulation, MipmapType mipmap, AnisotropicFilteringType anisotropicFiltering, TextureTransparencyType textureTransparency, std::vector<ObjMeshRender::Texture*> * texturePool, int updatePool)
{
  int textureMode = 0; // = OBJMESHRENDER_GL_MODULATE | OBJMESHRENDER_GL_NOMIPMAP | OBJMESHRENDER_GL_ANISOTROPICFILTERING

  switch(lightingModulation)
  {
    case REPLACE:
      textureMode |= OBJMESHRENDER_GL_REPLACE;
      break;
    case MODULATE:
      textureMode |= OBJMESHRENDER_GL_MODULATE;
      break;
  }

  switch(mipmap)
  {
    case USEMIPMAP:
      textureMode |= OBJMESHRENDER_GL_USEMIPMAP;
      break;
    case NOMIPMAP:
      textureMode |= OBJMESHRENDER_GL_NOMIPMAP;
      break;
  }

  switch(anisotropicFiltering)
  {
    case USEANISOTROPICFILTERING:
      textureMode |= OBJMESHRENDER_GL_USEANISOTROPICFILTERING;
      break;
    case NOANISOTROPICFILTERING:
      textureMode |= OBJMESHRENDER_GL_NOANISOTROPICFILTERING;
      break;
  }

  meshRender->loadTextures(textureMode, texturePool, updatePool);
  if (meshRender->numTextures() > 0)
    hasTextures_ = true;
  else
    hasTextures_ = false;

  EnableTextures();

  switch(textureTransparency)
  {
    case USETEXTURETRANSPARENCY:
      if (meshRender->maxBytesPerPixelInTextures() == 4)
        renderMode |= OBJMESHRENDER_TRANSPARENCY;
      break;     
    case NOTEXTURETRANSPARENCY:
      break;
  }

  return 0;
}

void SceneObject::RenderNormals()
{
  double normalLength = 0.1;
  meshRender->renderNormals(normalLength);
}

void SceneObject::BuildNeighboringStructure()
{
  mesh->buildVertexFaceNeighbors();
}

void SceneObject::ComputeMeshGeometricParameters(Vec3d * centroid, double * radius)
{
  mesh->getMeshGeometricParameters(centroid, radius);
}

void SceneObject::ComputeMeshRadius(Vec3d & centroid, double * radius)
{
  mesh->getMeshRadius(centroid, radius);
}

void SceneObject::ExportMeshGeometry(int * numVertices, double ** vertices, int * numTriangles, int ** triangles)
{
  mesh->exportGeometry(numVertices, vertices, numTriangles, triangles, NULL, NULL);
}

void SceneObject::PrintBitmapString(float x, float y, float z, const char* s)
{
  glRasterPos3f(x,y,z);
  if (s && strlen(s)) 
  {
    while (*s) 
    {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *s);
      s++;
    }
  }
}

void SceneObject::PrintBitmapInteger(float x, float y, float z, int i)
{
  char s[200];
  sprintf(s,"%d",i);
  PrintBitmapString(x,y,z,s);
}

void SceneObject::TransformRigidly(double * centerOfMass, double * R)
{
  Vec3d cv(centerOfMass);
  Mat3d Rv(R);
  mesh->transformRigidly(cv, Rv);
}

