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
#ifdef WIN32
  #include "GL/glew.h"
#endif
#include "objMeshGPUDeformer_coarseToFine.h"
#include "coarseToFine-shaders.cpp"

ObjMeshGPUDeformer_coarseToFine::ObjMeshGPUDeformer_coarseToFine() : ObjMeshGPUDeformer()
{
  interp_vertices = NULL;
  interp_weights = NULL;
  uCoarseTextureData = NULL;
}

void ObjMeshGPUDeformer_coarseToFine::Init(ObjMesh * mesh_, ObjMeshRender * meshRender_, int numCoarseVertices_, int interp_numElementVertices_, int * interp_vertices_, double * interp_weights_, int renderingMode_)
{
  ObjMeshGPUDeformer::Init(mesh_, meshRender_, renderingMode_);

  //internalFormat = GL_RGBA_FLOAT16_ATI;
  //internalFormat = GL_RGBA_FLOAT32_ATI;
  //internalFormat = GL_RGBA16F_ARB;
  internalFormat = GL_RGBA32F_ARB;

  numCoarseVertices = numCoarseVertices_;
  interp_numElementVertices = interp_numElementVertices_;

  if (interp_numElementVertices != 4)
  {
    printf("Error: ObjMeshGPUDeformer_coarseToFine only supports tet meshes.\n");
    exit(1);
  }

  interp_vertices = (int*) malloc (sizeof(int) * interp_numElementVertices * numVertices);
  memcpy(interp_vertices, interp_vertices_, sizeof(int) * interp_numElementVertices * numVertices);
  interp_weights = (double*) malloc (sizeof(double) * interp_numElementVertices * numVertices);
  memcpy(interp_weights, interp_weights_, sizeof(double) * interp_numElementVertices * numVertices);

  //sprintf(shaderFilename, "/Users/barbic/libraries/ObjMeshGPUDeformer_uUq/interp-shaders.cg");
  
  printf("Initializing ObjMeshGPUDeformer_coarseToFine...\n");
  printf ("Num vertices is %d.\n", numVertices);

  int code = InitializeCGShaders();
  if (code != 0)
  {
    //free(compilerOptions[0]);
    //free(compilerOptions);
    throw 32;
  }

  // ------- generate vertex texture coordinates ----------------------------
  InitVertexDeformationTexture();

  // ------- generate the coarse deformation texture ----------------------------
  InitCoarseDeformationTexture();

  // -------------- init the pbuffer -----------------------------------------

  if (InitRTT() != 0)
    throw 34;

  // ------- generate fragment U texture -------------------------------------

  InitInterpolationTexture();

  // load fragment programs
  EnableRTT();
  
  // fragment program init
  // the pass 1 fragment program that computes u=Uq
  
  /*
    Fragment_InterpolationProgram = cgCreateProgramFromFile(Context,
  	  CG_SOURCE, 
      //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_coarseToFine\\shaders.cg",
      shaderFilename,
  	  Fragment_InterpolationProfile,
  	  "Interpolation", NULL); // NULL for compiler options
  */
  
  Fragment_InterpolationProgram = cgCreateProgram(Context, CG_SOURCE, GPGPUInterpCode,
    Fragment_InterpolationProfile, "Interpolation", NULL); 
  
  if (Fragment_InterpolationProgram == NULL)
  {
    printf("Error loading Fragment_InterpolationProgram program.\n");
    throw 35;
  }
  
  cgGLLoadProgram(Fragment_InterpolationProgram);
  interpolationTextureParam = cgGetNamedParameter(Fragment_InterpolationProgram, "interpolationTexture");
  interpolationTextureSizeParam = cgGetNamedParameter(Fragment_InterpolationProgram, "interpolationTextureSize");
  coarseDeformationTextureParam = cgGetNamedParameter(Fragment_InterpolationProgram, "coarseDeformationTexture");
  coarseDeformationTextureSizeParam = cgGetNamedParameter(Fragment_InterpolationProgram, "coarseDeformationTextureSize");
  vertexDeformationTextureSizeParam = cgGetNamedParameter(Fragment_InterpolationProgram, "vertexDeformationTextureSize");
  
  DisableRTT();

  PrintGLerror("before making display lists");

  MakeDisplayLists(renderingMode);

  printf ("ObjMeshGPUDeformer_coarseToFine initialization complete.\n");
}

void ObjMeshGPUDeformer_coarseToFine::InitCoarseDeformationTexture()
{
  coarseDeformationTextureSize = 4;
  while (coarseDeformationTextureSize * coarseDeformationTextureSize < numCoarseVertices)
    coarseDeformationTextureSize *= 2; 

  // Create the dynamic texture (render-to-texture target)
  GLenum err;
  if ((err = glGetError()) != GL_NO_ERROR)
  {
    const GLubyte * errString = gluErrorString(err);
    printf("Warning at (1): %s \n",errString);
  }

  glGenTextures( 1, &coarseDeformationTextureID );
  glBindTexture( GL_TEXTURE_2D, coarseDeformationTextureID );

  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

  uCoarseTextureData = (float*) calloc (4 * coarseDeformationTextureSize * coarseDeformationTextureSize, sizeof(float));

  for(int i1 = 0 ; i1 < coarseDeformationTextureSize; i1++)
    for(int j1 = 0 ; j1 < coarseDeformationTextureSize; j1++)
    {
      uCoarseTextureData[(j1 * coarseDeformationTextureSize + i1) * 4 + 0] = 0.0;
      uCoarseTextureData[(j1 * coarseDeformationTextureSize + i1) * 4 + 1] = 2.0;
      uCoarseTextureData[(j1 * coarseDeformationTextureSize + i1) * 4 + 2] = 0.0;
      uCoarseTextureData[(j1 * coarseDeformationTextureSize + i1) * 4 + 3] = 1.0;
    }

  //free(test);
  //
  //if ((err = glGetError()) != GL_NO_ERROR)
  //{
    //const GLubyte * errString = gluErrorString(err);
    //printf("Warning at (2): %s \n",errString);
  //}

  glTexImage2D(GL_TEXTURE_2D, 0, internalFormat,
    coarseDeformationTextureSize, coarseDeformationTextureSize, 0, 
    GL_RGBA, GL_FLOAT, uCoarseTextureData);

  if ((err = glGetError()) != GL_NO_ERROR)
  {
    const GLubyte * errString = gluErrorString(err);
    printf("Warning at (3): %s \n", errString);
  }

  printf("Coarse deformation texture initialized. Size of texture is %d x %d .\n", coarseDeformationTextureSize, coarseDeformationTextureSize);
}

void ObjMeshGPUDeformer_coarseToFine::InitInterpolationTexture()
{
  interpolationTextureSize = 4;
  while (interpolationTextureSize * interpolationTextureSize < 2 * numVertices) // 4 floats per pixel
    interpolationTextureSize *= 2; 

  float * interpolationTextureData = (float*) calloc (4 * interpolationTextureSize * interpolationTextureSize, sizeof(float));

  int row = 0;
  int column = 0;
  for(int i=0; i<numVertices; i++)
  {
    for(int j=0; j<4; j++)
    {
      interpolationTextureData[4 * (row * interpolationTextureSize + column) + j] = (float)(interp_vertices[4*i+j]);
      //interpolationTextureData[4 * (row * interpolationTextureSize + column) + j] = 5.0; 
    }

    for(int j=0; j<4; j++)
    {
      interpolationTextureData[4* (row * interpolationTextureSize + column + 1) + j] = interp_weights[4*i+j];
      //interpolationTextureData[4* (row * interpolationTextureSize + column + 1) + j] = -5.0; 
    }

    column += 2;
    if (column >= interpolationTextureSize)
    {
      column = 0;
      row++;
    }
  }

  glGenTextures(1, &interpolationTextureID);
  glBindTexture(GL_TEXTURE_2D, interpolationTextureID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // GL_NEAREST_MIPMAP_NEAREST
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glTexImage2D(GL_TEXTURE_2D, 0, internalFormat,
      interpolationTextureSize, interpolationTextureSize, 0,
      GL_RGBA, GL_FLOAT, interpolationTextureData);

  invInterpolationTextureSize = 1.0 / interpolationTextureSize;
  free(interpolationTextureData);

  printf("Interpolation texture initialized. Size is %d x %d .\n", interpolationTextureSize, interpolationTextureSize);
}

ObjMeshGPUDeformer_coarseToFine::~ObjMeshGPUDeformer_coarseToFine()
{
  DeleteCGShaders();

  if (Fragment_InterpolationProgram)
    cgDestroyProgram(Fragment_InterpolationProgram);

  if (interp_vertices != NULL)
    free(interp_vertices);

  if (interp_weights != NULL)
    free(interp_weights);

  if (uCoarseTextureData != NULL)
    free(uCoarseTextureData);

  DeleteRTT();
}

int ObjMeshGPUDeformer_coarseToFine::InitializeCGShaders()
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

  cgSetErrorCallback(ObjMeshGPUDeformer_coarseToFine::cgErrorCallback);
  Context = cgCreateContext();

  VertexProfile = cgGLGetLatestProfile(CG_GL_VERTEX);
  cgGLSetOptimalOptions(VertexProfile);
  FragmentProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(FragmentProfile);
  Fragment_InterpolationProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(Fragment_InterpolationProfile);

  //printf("Loading shader from %s.\n", shaderFilename);


    // add pass 2 vertex shader
/*
  VertexPass2Program = cgCreateProgramFromFile(Context,
    CG_SOURCE, //"D:/barbic/libraries/ObjMeshGPUDeformer_coarseToFine/shaders.cg",
      shaderFilename, VertexProfile, "vertexShaderPass2", NULL);
*/

  VertexPass2Program = cgCreateProgram(Context, CG_SOURCE, GPGPUInterpCode, VertexProfile, "vertexShaderPass2", NULL);

  if(VertexPass2Program != NULL)
  {
    // Vertex shader only needs to be loaded once 
    cgGLLoadProgram(VertexPass2Program);
    // Bind parameters to give access to variables in the shader 
    //texDeltaParam = cgGetNamedParameter(VertexPass2Program, "texDelta");
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

  /*
    VertexPass2ProgramShadow = cgCreateProgramFromFile(Context,
      CG_SOURCE, //"D:/barbic/libraries/ObjMeshGPUDeformer_coarseToFine/shaders.cg",
      shaderFilename, VertexProfile, "vertexShaderShadowPass2", NULL);
  */
    
  VertexPass2ProgramShadow = cgCreateProgram(Context,
    CG_SOURCE, GPGPUInterpCode, VertexProfile, "vertexShaderShadowPass2", NULL);

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

  // create the vertex program for points
/*
  VertexPass2ProgramPoints = cgCreateProgramFromFile(Context,
  	CG_SOURCE, //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_coarseToFine\\shaders.cg",
    "shaders.cg", VertexProfile, "vertexShader_Points_Pass2", (const char **)compilerOptions);
    */

/*
  VertexPass2ProgramPoints = cgCreateProgram(Context,
    CG_SOURCE, GPGPUInterpCode, VertexProfile,
    "vertexShader_Points_Pass2", (const char **)compilerOptions);

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
*/

  // create the vertex program for rendering edges
  //VertexPass2ProgramEdges = cgCreateProgramFromFile(Context,
    //CG_SOURCE, 
    ////"D:\\barbic\\libraries\\ObjMeshGPUDeformer_coarseToFine\\shaders.cg",
    //"shaders.cg",
    //VertexProfile,
    //"vertexShader_Edges_Pass2", (const char **)compilerOptions);

/*
  VertexPass2ProgramEdges = cgCreateProgram(Context,
    CG_SOURCE, GPGPUInterpCode,
    VertexProfile,
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
*/

/*
  // add pass 2 fragment shader 
  FragmentPass2Program = cgCreateProgramFromFile(Context,
  	CG_SOURCE, //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_coarseToFine\\shaders.cg",
    "interp_shaders.cg", FragmentProfile, "fragmentShaderPass2", NULL);
*/

  // add pass 2 fragment shader 
  FragmentPass2Program = cgCreateProgram(Context,
    CG_SOURCE, GPGPUInterpCode, FragmentProfile, "fragmentShaderPass2", NULL);

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

/*
  // add pass 2 no texture fragment shader 
  FragmentPass2Program = cgCreateProgramFromFile(Context,
  	CG_SOURCE, 
    //"D:\\barbic\\libraries\\ObjMeshGPUDeformer_coarseToFine\\shaders.cg",
    "interp_shaders.cg", FragmentProfile, "fragmentShaderPass2NoTexture", NULL);
*/

  // add pass 2 fragment shader 
  FragmentPass2ProgramNoTexture = cgCreateProgram(Context,
    CG_SOURCE, GPGPUInterpCode, FragmentProfile, "fragmentShaderPass2NoTexture", NULL);

  if(FragmentPass2ProgramNoTexture != NULL)
  {
    // Fragment shader only needs to be loaded once 
    cgGLLoadProgram(FragmentPass2ProgramNoTexture);
  }
  else
  {
    printf("Error loading pass 2 no texture fragment program.\n");
    return 1;
  }

  hasStaticBeenInitialized = true;

  return 0;
}

void ObjMeshGPUDeformer_coarseToFine::SetCoarseDeformations(double * uCoarse)
{
  for(int i=0; i<numCoarseVertices; i++)
  {
    for(int j=0; j<3; j++)
    {
      uCoarseTextureData[4 * i + j] = (float)(uCoarse[3*i+j]);
      //printf("%G ", uCoarseTextureData[4 * i + j]);
    }
    //printf("\n");
    uCoarseTextureData[4 * i + 3] = 0.0;
  }

  UploadDeformations();

  //PrintCoarseDeformationTexture();

  InterpolateDeformations(); 

  //PrintFineDeformationTexture();
  //exit(1); 
}

void ObjMeshGPUDeformer_coarseToFine::UploadDeformations()
{
  glBindTexture(GL_TEXTURE_2D, coarseDeformationTextureID);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, coarseDeformationTextureSize, coarseDeformationTextureSize, GL_RGBA, GL_FLOAT, uCoarseTextureData);

  //glTexImage2D( GL_TEXTURE_2D, 0, internalFormat,
    //coarseDeformationTextureSize, coarseDeformationTextureSize, 0, 
    //GL_RGBA, GL_FLOAT, uCoarseTextureData );
}

void ObjMeshGPUDeformer_coarseToFine::InterpolateDeformations()
{
  EnableRTT();

  cgGLEnableProfile(Fragment_InterpolationProfile);
  cgGLBindProgram(Fragment_InterpolationProgram);
  
  cgGLSetTextureParameter(interpolationTextureParam, interpolationTextureID);
  cgGLEnableTextureParameter(interpolationTextureParam);

  cgGLSetTextureParameter(coarseDeformationTextureParam, coarseDeformationTextureID);
  cgGLEnableTextureParameter(coarseDeformationTextureParam);
 
  cgGLSetParameter1f(vertexDeformationTextureSizeParam, vertexDeformationTextureSize);
  cgGLSetParameter1f(coarseDeformationTextureSizeParam, coarseDeformationTextureSize);
  cgGLSetParameter1f(interpolationTextureSizeParam, interpolationTextureSize);

  float maxX = 1.0;
  // don't draw the upper image part where there is no data
  int maxRow = (numVertices-1) / vertexDeformationTextureSize; // maxRow can take values 0...(vertexDeformationTextureSize-1)
  // maxY is upper side of row maxRow
  float maxY = 1.0 * (maxRow+1) / vertexDeformationTextureSize;
  //printf("maxY=%G\n", maxY);
  //maxY = 1.0; 

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

  //glActiveTextureARB(GL_TEXTURE0_ARB); 

  // render the quad, size is vertexDeformationTextureSize x vertexDeformationTextureSize
  //glDisable(GL_LIGHTING);
  //glColor3f(1,0,1);
  glBegin(GL_QUADS);
    glMultiTexCoord2fARB(GL_TEXTURE1_ARB, 0, 0);
    glVertex3f(0,0,0);

    glMultiTexCoord2fARB(GL_TEXTURE1_ARB, maxX, 0);
    //glTexCoord2f(maxX,0); 
    glVertex3f(maxX,0,0);

    glMultiTexCoord2fARB(GL_TEXTURE1_ARB, maxX, maxY);
    //glTexCoord2f(maxX,maxY); 
    glVertex3f(maxX,maxY,0);

    glMultiTexCoord2fARB(GL_TEXTURE1_ARB, 0, maxY);
    //glTexCoord2f(0,maxY);
    glVertex3f(0,maxY,0);
  glEnd();

  cgGLDisableTextureParameter(interpolationTextureParam);
  cgGLDisableTextureParameter(coarseDeformationTextureParam);

  cgGLDisableProfile(Fragment_InterpolationProfile);

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);

  DisableRTT();
}

void ObjMeshGPUDeformer_coarseToFine::PrintFineDeformationTexture()
{
  PrintGLerror("Entering readback");
  float * buffer = (float*) malloc (sizeof(float) * 4 * vertexDeformationTextureSize * vertexDeformationTextureSize);

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, buffer);

  for(int i=0; i<numVertices; i++)
    printf("%G %G %G\n", buffer[4*i+0], buffer[4*i+1], buffer[4*i+2]);
  free(buffer);
  PrintGLerror("Leaving readback");
}

void ObjMeshGPUDeformer_coarseToFine::PrintCoarseDeformationTexture()
{
  PrintGLerror("Entering readback");
  float * buffer = (float*) malloc (sizeof(float) * 4 * coarseDeformationTextureSize * coarseDeformationTextureSize);

  glBindTexture(GL_TEXTURE_2D, coarseDeformationTextureID);
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, buffer);

  printf("Coarse texture:\n");
  for(int i=0; i<3; i++)
    printf("%G %G %G\n", buffer[4*i+0], buffer[4*i+1], buffer[4*i+2]);
  free(buffer);
  PrintGLerror("Leaving readback");
}

// initialize static members
bool ObjMeshGPUDeformer_coarseToFine::hasStaticBeenInitialized = false;

CGprofile ObjMeshGPUDeformer_coarseToFine::Fragment_InterpolationProfile;
CGprogram ObjMeshGPUDeformer_coarseToFine::Fragment_InterpolationProgram;
CGprogram ObjMeshGPUDeformer_coarseToFine::FragmentPass2ProgramNoTexture;

