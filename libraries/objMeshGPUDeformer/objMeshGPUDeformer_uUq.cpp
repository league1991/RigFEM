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
#include "stdafx.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "objMeshGPUDeformer_uUq.h"
#include "uUq-shaders.cpp"

extern char fragment_uUqShaderProgramCode [];
extern char fragmentShaderPass2ProgramCode [];
extern char vertexShaderProgramCode [];

ObjMeshGPUDeformer_uUq::ObjMeshGPUDeformer_uUq() : ObjMeshGPUDeformer() {}

void ObjMeshGPUDeformer_uUq::Init(ObjMesh * mesh_, ObjMeshRender * meshRender_, int r_, double * U_, int renderingMode_)
{
  ObjMeshGPUDeformer::Init(mesh_, meshRender_, renderingMode_);

  r = r_;
  U = U_;
  q = (float*) calloc (sizeof(float), r);
  printf("Initializing ObjMeshGPUDeformer_uUq...\n");
  printf ("r is %d, num vertices is %d\n", r, numVertices);

  // init compiler options
  compilerOptions = (char**) malloc (sizeof(char*) * 2);
  compilerOptions[0] = (char*) malloc (sizeof(char) * 20);
  sprintf(compilerOptions[0], "-DR=%d%c", r, 0);
  compilerOptions[1] = NULL;

  int code = InitializeCGShaders();
  if (code != 0)
  {
    free(compilerOptions[0]);
    free(compilerOptions);
    throw 32;
  }

  // -------------- init the pbuffer -----------------------------------------

  if (InitRTT() != 0)
    throw 34;

  // ------- generate fragment U texture -------------------------------------

  EnableRTT();
    InitUTexture();
  DisableRTT();

  // load fragment programs

  Fragment_uUqProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(Fragment_uUqProfile);

  EnableRTT();

  // fragment program init
  // the pass 1 fragment program that computes u=Uq

  printf("Init pass 1 fragment program.\n");

  /*
  Fragment_uUqProgram = cgCreateProgramFromFile(Context,
  	CG_SOURCE, 
    //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_uUq\\shaders.cg",
    "shaders.cg",
  	Fragment_uUqProfile,
  	"uUq", (const char **)compilerOptions); // used to be NULL for options
  */

  Fragment_uUqProgram = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode, Fragment_uUqProfile,
    "uUq", (const char **)compilerOptions); // used to be NULL for options

  if (Fragment_uUqProgram == NULL)
  {
    printf("Error loading Fragment_uUqProgram program.\n");
    throw 35;
  }

  printf("Pass 1 fragment program initialized.\n");

  cgGLLoadProgram(Fragment_uUqProgram);
  UtexParam = cgGetNamedParameter(Fragment_uUqProgram, "Utex");
  qParam = cgGetNamedParameter(Fragment_uUqProgram, "q");
  vertexDeformationTextureSizeParam = cgGetNamedParameter(Fragment_uUqProgram, "size");
  SIZEParam = cgGetNamedParameter(Fragment_uUqProgram, "SIZE");
  SIZEdivRParam = cgGetNamedParameter(Fragment_uUqProgram, "SIZEdivR");
  texUDeltaParam = cgGetNamedParameter(Fragment_uUqProgram, "texUDelta");

  DisableRTT();

  PrintGLerror("before making display lists");

  MakeDisplayLists(renderingMode);

  printf ("ObjMeshGPUDeformer_uUq initialization complete.\n");
}

ObjMeshGPUDeformer_uUq::~ObjMeshGPUDeformer_uUq()
{
  DeleteCGShaders();
  if (Fragment_uUqProgram)
    cgDestroyProgram(Fragment_uUqProgram);
  DeleteRTT();
}

void ObjMeshGPUDeformer_uUq::InitUTexture()
{
  UTextureSize = 1 << (1+(int)(0.5*log(1.0f*r*numVertices)/log (2.0f)));
  if (r*numVertices > UTextureSize * ((UTextureSize / r) * r))
    UTextureSize *= 2;

  float * fragmentUTextureData = (float*) calloc (4 * UTextureSize * UTextureSize, sizeof(float));

  int row = 0;
  int column = 0;
  for(int i=0; i< numVertices; i++)
  {
    if (column + r >= UTextureSize)
    {
      column = 0;
      row++;
    }

    for(int j=0; j<r; j++)
    {
      fragmentUTextureData[4 * (row * UTextureSize + column) + 0] = U[3*numVertices*j + 3*i+0];
      fragmentUTextureData[4 * (row * UTextureSize + column) + 1] = U[3*numVertices*j + 3*i+1];
      fragmentUTextureData[4 * (row * UTextureSize + column) + 2] = U[3*numVertices*j + 3*i+2];
      fragmentUTextureData[4 * (row * UTextureSize + column) + 3] = 1.0;
      column++;
    }
  }

  glGenTextures(1, &UTextureID);
  glBindTexture(GL_TEXTURE_2D, UTextureID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // GL_NEAREST_MIPMAP_NEAREST

  //int internalFormat = GL_RGBA_FLOAT16_ATI;
  //int internalFormat = GL_RGBA_FLOAT32_ATI;
  int internalFormat = GL_RGBA16F_ARB;
  //int internalFormat = GL_RGBA32F_ARB;
  glTexImage2D(GL_TEXTURE_2D, 0, internalFormat,
      UTextureSize, UTextureSize, 0,
      GL_RGBA, GL_FLOAT, fragmentUTextureData);

  invFragmentTextureSize = 1.0 / UTextureSize;
  numPixelsInARow = UTextureSize / r;
  free(fragmentUTextureData);

  printf("Modal matrix texture initialized. Size is %d x %d .\n", UTextureSize, UTextureSize);
}

int ObjMeshGPUDeformer_uUq::InitializeCGShaders()
{
  if (hasStaticBeenInitialized)
    return 0;

  // -------------- init Cg ---------------------------------------------------

  // check for presence of vertex texture fetches
  if (!glh_extension_supported("GL_NV_vertex_program3"))
  {
    printf("Error: did not detect the vertex texture fetch NV_vertex_program3 extension.\n");
    return 1;
  }

  cgSetErrorCallback(ObjMeshGPUDeformer_uUq::cgErrorCallback);
  Context = cgCreateContext();

  VertexProfile = cgGLGetLatestProfile(CG_GL_VERTEX);
  cgGLSetOptimalOptions(VertexProfile);
  FragmentProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(FragmentProfile);

  printf("Initializing pass 2 vertex shader.\n");
  // add pass 2 vertex shader
/*
  VertexProgram = cgCreateProgramFromFile(Context,
    CG_SOURCE, 
    //"D:/barbic/libraries/ObjMeshGPUDeformer_uUq/shaders.cg",
    "shaders.cg",
  	VertexProfile,
  	"vertexShaderPass2", (const char **)compilerOptions);
*/
  VertexPass2Program = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode, VertexProfile,
    "vertexShaderPass2", (const char **)compilerOptions);

  if(VertexPass2Program != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2Program);
    // Bind parameters to give access to variables in the shader 
    //texParam = cgGetNamedParameter(VertexProgram, "tex");
    //texDeltaParam = cgGetNamedParameter(VertexProgram, "texDelta");
    KaParam = cgGetNamedParameter(VertexPass2Program, "Ka");
    KdParam = cgGetNamedParameter(VertexPass2Program, "Kd");
    KsParam = cgGetNamedParameter(VertexPass2Program, "Ks");
    shininessParam = cgGetNamedParameter(VertexPass2Program, "shininess");
    ModelViewProjParam = cgGetNamedParameter(VertexPass2Program, "ModelViewProj");
    ModelViewITParam = cgGetNamedParameter(VertexPass2Program, "ModelViewIT");
    LightPos1Param = cgGetNamedParameter(VertexPass2Program, "LightPos1");
    LightPos2Param = cgGetNamedParameter(VertexPass2Program, "LightPos2");
    LightPos3Param = cgGetNamedParameter(VertexPass2Program, "LightPos3");
    LightPos4Param = cgGetNamedParameter(VertexPass2Program, "LightPos4");
    Light1IntensityParam = cgGetNamedParameter(VertexPass2Program, "Light1Intensity");
    Light2IntensityParam = cgGetNamedParameter(VertexPass2Program, "Light2Intensity");
    Light3IntensityParam = cgGetNamedParameter(VertexPass2Program, "Light3Intensity");
    Light4IntensityParam = cgGetNamedParameter(VertexPass2Program, "Light4Intensity");
    AmbientIntensityParam = cgGetNamedParameter(VertexPass2Program, "AmbientIntensity");
    //vertexTextureCoordParam = cgGetNamedParameter(VertexPass2Program, "IN.vertexTextureCoord");
  }
  else
  {
    printf("Error loading vertex program.\n");
    return 1;
  }

  printf("Initializing pass 2 vertex shadow shader.\n");

  VertexPass2ProgramShadow = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode,
    VertexProfile, "vertexShaderShadowPass2", (const char **)compilerOptions);

  if(VertexPass2ProgramShadow != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2ProgramShadow);
    // Bind parameters to give access to variables in the shader 
    ShadowIntensityParam = cgGetNamedParameter(VertexPass2ProgramShadow, "ShadowIntensity");
    ModelViewProjShadowParam = cgGetNamedParameter(VertexPass2ProgramShadow, "ModelViewProj");
  }
  else
  {
    printf("Error loading vertex shadow program.\n");
    return 1;
  }
/*
  // load pass 2 vertex shader with normals
  VertexPass2ProgramDeformedNormals = cgCreateProgram(Context,
    CG_SOURCE, GPGPUuUqCode,
    VertexProfile,
    "vertexShaderPass2WithDefoNormals", (const char **)compilerOptions);

  if(VertexPass2ProgramDeformedNormals != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2ProgramDeformedNormals);
    // Bind parameters to give access to variables in the shader 
    //texParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "tex");
    //texDeltaParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "texDelta");
    KaDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Ka");
    KdDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Kd");
    KsDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Ks");
    shininessDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "shininess");
    ModelViewProjDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "ModelViewProj");
    ModelViewITDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "ModelViewIT");
    LightPos1DeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "LightPos1");
    LightPos2DeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "LightPos2");
    LightPos3DeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "LightPos3");
    LightPos4DeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "LightPos4");
    Light1IntensityDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Light1Intensity");
    Light2IntensityDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Light2Intensity");
    Light3IntensityDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Light3Intensity");
    Light4IntensityDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "Light4Intensity");
    AmbientIntensityDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "AmbientIntensity");
    //vertexTextureCoordDeformedNormalsParam = cgGetNamedParameter(VertexPass2ProgramDeformedNormals, "IN.vertexTextureCoord");
  }
  else
  {
    printf("Error loading vertex program.\n");
    return 1;
  }
*/

  // create the vertex program for points
/*
  VertexPass2ProgramPoints = cgCreateProgramFromFile(Context, CG_SOURCE, 
    //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_uUq\\shaders.cg",
    "shaders.cg", VertexProfile,
    "vertexShader_Points_Pass2", (const char **)compilerOptions);
*/

  printf("Initializing pass 2 vertex point shader.\n");

  VertexPass2ProgramPoints = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode,
    VertexProfile, "vertexShader_Points_Pass2", (const char **)compilerOptions);

  if(VertexPass2ProgramPoints != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2ProgramPoints);
    // Bind parameters to give access to variables in the shader 
    //texPointParam = cgGetNamedParameter(VertexPass2ProgramPoints, "tex");
    ModelViewProjPointParam = cgGetNamedParameter(VertexPass2ProgramPoints, "ModelViewProj");
    ModelViewITPointParam = cgGetNamedParameter(VertexPass2ProgramPoints, "ModelViewIT");
  }
  else
  {
    printf("Error loading vertex program for points.\n");
    return 1;
  }

  // create the vertex program for rendering edges
/*
  VertexPass2ProgramEdges = cgCreateProgramFromFile(Context, CG_SOURCE, 
    //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_uUq\\shaders.cg",
    "shaders.cg", VertexProfile,
    "vertexShader_Edges_Pass2", (const char **)compilerOptions);
*/

  printf("Initializing pass 2 vertex edge shader.\n");

  VertexPass2ProgramEdges = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode, VertexProfile,
    "vertexShader_Edges_Pass2", (const char **)compilerOptions);

  if(VertexPass2ProgramEdges != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2ProgramEdges);
    // Bind parameters to give access to variables in the shader 
    //texEdgeParam = cgGetNamedParameter(VertexPass2ProgramEdges, "tex");
    ModelViewProjEdgeParam = cgGetNamedParameter(VertexPass2ProgramEdges, "ModelViewProj");
    ModelViewITEdgeParam = cgGetNamedParameter(VertexPass2ProgramEdges, "ModelViewIT");
  }
  else
  {
    printf("Error loading vertex program for edges.\n");
    return 1;
  }

/*
  // add pass 2 fragment shader 
  FragmentPass2Program = cgCreateProgramFromFile(Context, CG_SOURCE, 
    //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_uUq\\shaders.cg",
    "shaders.cg", FragmentProfile,
    "fragmentShaderPass2", (const char **)compilerOptions);
*/

  printf("Initializing pass 2 fragment shader.\n");

  // add pass 2 fragment shader 
  FragmentPass2Program = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode, FragmentProfile,
    "fragmentShaderPass2", (const char **)compilerOptions);

  if(FragmentPass2Program != NULL)
  {
    // Fragment shader only needs to be loaded once 
    cgGLLoadProgram(FragmentPass2Program);
    // Bind parameters to give access to variables in the shader 
    texPass2Param = cgGetNamedParameter(FragmentPass2Program, "texPass2");
  }
  else
  {
    printf("Error loading pass 2 fragment program.\n");
    return 1;
  }

  printf("end of static init.\n");

  hasStaticBeenInitialized = true;

  return 0;
}

void ObjMeshGPUDeformer_uUq::Setqfv(float * q_)
{
  memcpy(q,q_,sizeof(float)*r);
}

void ObjMeshGPUDeformer_uUq::Setqdv(double * q_)
{
  for(int i=0; i<r; i++)
    q[i] = (float)q_[i];
}

void ObjMeshGPUDeformer_uUq::Compute_uUq()
{
  EnableRTT();

  cgGLEnableProfile(Fragment_uUqProfile);
  cgGLBindProgram(Fragment_uUqProgram);

  cgGLSetTextureParameter(UtexParam, UTextureID);
  cgGLSetParameter1f(vertexDeformationTextureSizeParam, vertexDeformationTextureSize);
  cgGLSetParameter1f(SIZEParam, UTextureSize);
  cgGLSetParameter1f(texUDeltaParam,invFragmentTextureSize);
  cgGLSetParameter1f(SIZEdivRParam, numPixelsInARow);
  cgGLEnableTextureParameter(UtexParam);

  cgGLSetParameterArray1f(qParam, 0, r, q);

  float maxX = 1.0;
  // don't draw the upper image part where there is no data
  int maxRow = (numVertices-1) / vertexDeformationTextureSize; // maxRow can take values 0...(vertexDeformationTextureSize-1)
  // maxY is upper side of row maxRow
  float maxY = 1.0 * (maxRow+1) / vertexDeformationTextureSize;

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, 1, 0, 1, 0, 1);
  glViewport(0, 0, vertexDeformationTextureSize, vertexDeformationTextureSize);

  // render the quad, size is vertexDeformationTextureSize x vertexDeformationTextureSize
  //glDisable(GL_LIGHTING);
  //glColor3f(1,0,1);
  glBegin(GL_QUADS);
    glTexCoord2f(0,0);
    glVertex3f(0,0,0);

    glTexCoord2f(maxX,0);
    glVertex3f(maxX,0,0);

    glTexCoord2f(maxX,maxY);
    glVertex3f(maxX,maxY,0);

    glTexCoord2f(0,maxY);
    glVertex3f(0,maxY,0);
  glEnd();

  cgGLDisableTextureParameter(UtexParam);
  cgGLDisableProfile(Fragment_uUqProfile);

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

  DisableRTT();

  glFinish();
}

void ObjMeshGPUDeformer_uUq::Clone(ObjMeshGPUDeformer_uUq * ObjMeshGPUDeformer_uUq_source)
{
  mesh = ObjMeshGPUDeformer_uUq_source->mesh;
  numVertices = ObjMeshGPUDeformer_uUq_source->numVertices;
  r = ObjMeshGPUDeformer_uUq_source->r;
  numTriangles = ObjMeshGPUDeformer_uUq_source->numTriangles;
  numGroups = ObjMeshGPUDeformer_uUq_source->numGroups;
  U = ObjMeshGPUDeformer_uUq_source->U;

  numGroupTriangles = (int*) malloc (sizeof(int) * numGroups);
  for(int i=0; i<numGroups; i++)
    numGroupTriangles[i] = (ObjMeshGPUDeformer_uUq_source->numGroupTriangles)[i];

  compilerOptions = (char**) malloc (sizeof(char*) * 2);
  compilerOptions[0] = (char*) malloc (sizeof(char) * 20);
  sprintf(compilerOptions[0],"-DR=%d%c",r,0);
  compilerOptions[1] = NULL;

  q = (float*) calloc (r, sizeof(float));
  memcpy(q, ObjMeshGPUDeformer_uUq_source->q, sizeof(float) * r);

  displayListStart = ObjMeshGPUDeformer_uUq_source->displayListStart;
  displayListPoints = ObjMeshGPUDeformer_uUq_source->displayListPoints;
  displayListEdgesStart = ObjMeshGPUDeformer_uUq_source->displayListEdgesStart;
  numGroups = ObjMeshGPUDeformer_uUq_source->numGroups;

  UTextureID = ObjMeshGPUDeformer_uUq_source->UTextureID;
  renderingMode = ObjMeshGPUDeformer_uUq_source->renderingMode;

  memcpy(lightPos, ObjMeshGPUDeformer_uUq_source->lightPos, sizeof(float) * 16);
  memcpy(lightIntensity, ObjMeshGPUDeformer_uUq_source->lightIntensity, sizeof(float) * 4);
  ambientIntensity = ObjMeshGPUDeformer_uUq_source->ambientIntensity;
  ambientEnabled = ObjMeshGPUDeformer_uUq_source->ambientEnabled;
  diffuseEnabled = ObjMeshGPUDeformer_uUq_source->diffuseEnabled;
  specularEnabled = ObjMeshGPUDeformer_uUq_source->specularEnabled;

  UTextureSize = ObjMeshGPUDeformer_uUq_source->UTextureSize;
  invFragmentTextureSize = ObjMeshGPUDeformer_uUq_source->invFragmentTextureSize;
  numPixelsInARow = ObjMeshGPUDeformer_uUq_source->numPixelsInARow;

  InitVertexDeformationTexture();
  free(gpgpuVertexTextureCoordinates);
  SetDerivedData(ObjMeshGPUDeformer_uUq_source->GetDerivedData());

  Fragment_uUqProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(Fragment_uUqProfile);

  EnableRTT();

  Fragment_uUqProgram = cgCreateProgram(Context, CG_SOURCE, GPGPUuUqCode,
    Fragment_uUqProfile, "uUq", (const char **)compilerOptions); // used to be NULL for options

  if (Fragment_uUqProgram == NULL)
  {
    printf("Error loading Fragment_uUqProgram program.\n");
    throw 35;
  }

  DisableRTT();  

  cgGLLoadProgram(Fragment_uUqProgram);
  UtexParam = cgGetNamedParameter(Fragment_uUqProgram, "Utex");
  qParam = cgGetNamedParameter(Fragment_uUqProgram, "q");
  vertexDeformationTextureSizeParam = cgGetNamedParameter(Fragment_uUqProgram, "size");
  SIZEParam = cgGetNamedParameter(Fragment_uUqProgram, "SIZE");
  SIZEdivRParam = cgGetNamedParameter(Fragment_uUqProgram, "SIZEdivR");
  texUDeltaParam = cgGetNamedParameter(Fragment_uUqProgram, "texUDelta");

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    vboID = ObjMeshGPUDeformer_uUq_source->vboID;
    vboEdgesID = ObjMeshGPUDeformer_uUq_source->vboEdgesID;
  #endif

  vboNormalEnabled = ObjMeshGPUDeformer_uUq_source->vboNormalEnabled;
  vboTex1Enabled = ObjMeshGPUDeformer_uUq_source->vboTex1Enabled;
}

bool ObjMeshGPUDeformer_uUq::hasStaticBeenInitialized = false;

