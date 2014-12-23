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

#ifndef _OBJMESHGPUDEFORMER__H_
#define _OBJMESHGPUDEFORMER__H_

#ifdef WIN32
  #include <windows.h>
//#include <gl/glew.h>
#else
  #define GL_GLEXT_PROTOTYPES 1
#endif

#include "openGL-headers.h"

#if defined(WIN32) || defined(linux)
  #include <GL/glext.h>
#elif defined(__APPLE__)
  #include <OpenGL/glext.h>
#endif

#include <Cg/cg.h>
#include <Cg/cgGL.h>

#include "objMesh.h"
#include "objMeshRender.h"

#define OBJMESHGPUDEFORMER_USING_VBOS

// abstract virtual class, do not initialize directly
class ObjMeshGPUDeformer
{
public:
  ObjMeshGPUDeformer();
  void Init(ObjMesh * mesh, ObjMeshRender * meshRender, int renderingMode);
  virtual ~ObjMeshGPUDeformer() = 0;

  void Render();
  void RenderVertices();
  void RenderEdges();
  void RenderShadow(double intensity);

  void ReadBack_u(float * u);
  void ReadBack_u(double * u);

  void SetLightPosition(int lightID, float pos[4]); // note: fourth component will be ignored and always set to 1
  void SetLightIntensity(int lightID, float intensity);
  void SetAmbientIntensity(float intensity);

  void EnableAmbientLight(int ambientEnabled) { this->ambientEnabled = ambientEnabled; }
  void EnableDiffuseLight(int diffuseEnabled) { this->diffuseEnabled = diffuseEnabled; }
  void EnableSpecularLight(int specularEnabled) { this->specularEnabled = specularEnabled; }

protected:

  ObjMesh * mesh;
  ObjMeshRender * meshRender;
  int numVertices;
  int numGroups;
  int numTriangles;
  int * numGroupTriangles;

  int renderingMode;
  int displayListStart;
  int displayListEdgesStart;
  int displayListPoints;

  float lightPos[16]; // supports up to 4 lights
  float lightIntensity[4]; // supports up to 4 lights
  float ambientIntensity;
  bool ambientEnabled;
  bool diffuseEnabled;
  bool specularEnabled;

  float * gpgpuVertexTextureCoordinates;
  GLuint vertexDeformationTextureID;
  int vertexDeformationTextureSize;
  void InitVertexDeformationTexture();

  void MakeDisplayLists(int mode);
  void MakeDisplayListsTriangles(int mode);
  void MakeDisplayListsPoints();
  void MakeDisplayListsEdges();

  void DeleteCGShaders();

  // mode = 0: Render
  // mode = 1: RenderEdges
  void RenderMaster(int masterMode, void * data = NULL);

  static int glh_extension_supported(const char * extension);
  static void cgErrorCallback(void);
  void PrintGLerror(const char * msg);

  virtual void EnableRTT() = 0;
  virtual void DisableRTT() = 0;
  virtual void BindRT() = 0;
  virtual void UnbindRT() = 0;
  virtual int InitRTT() = 0;
  virtual void DeleteRTT() {}

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    GLuint * vboID;
    GLuint * vboEdgesID;
  #endif

  int * vboNormalEnabled;
  int * vboTex1Enabled;

  #ifdef WIN32
    void heap_check_();
  #endif

  static CGcontext Context; // needed to be made static for the purposes of the error callback function

  static CGprofile VertexProfile;// = CG_PROFILE_VP40;
  static CGprogram VertexPass2Program;
  static CGprogram VertexPass2ProgramShadow;
  //static CGprogram VertexPass2ProgramDeformedNormals;
  static CGprogram VertexPass2ProgramPoints;
  static CGprogram VertexPass2ProgramEdges;

  static CGprofile FragmentProfile; // = CG_PROFILE_FP30;
  static CGprogram FragmentPass2Program;

  // vertex shader parameters
  static CGparameter texParam;
  static CGparameter fragmentTexParam;
  static CGparameter KaParam;
  static CGparameter KdParam;
  static CGparameter KsParam;
  static CGparameter shininessParam;
  static CGparameter ModelViewProjParam;
  static CGparameter ModelViewITParam;
  static CGparameter LightPos1Param;
  static CGparameter LightPos2Param;
  static CGparameter LightPos3Param;
  static CGparameter LightPos4Param;
  static CGparameter Light1IntensityParam;
  static CGparameter Light2IntensityParam;
  static CGparameter Light3IntensityParam;
  static CGparameter Light4IntensityParam;
  static CGparameter AmbientIntensityParam;
  static CGparameter vertexTextureCoordParam;

  // vertex shader shadow
  static CGparameter ModelViewProjShadowParam;
  static CGparameter ShadowIntensityParam;

  // vertex shader for points parameters
  static CGparameter ModelViewProjPointParam;
  static CGparameter ModelViewITPointParam;

  // vertex shader for edges parameters
  static CGparameter ModelViewProjEdgeParam;
  static CGparameter ModelViewITEdgeParam;

  // pass 2 fragment shader parameters
  static CGparameter texPass2Param;
};

#endif

