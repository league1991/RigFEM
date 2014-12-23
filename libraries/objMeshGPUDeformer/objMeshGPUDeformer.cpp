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
#include "objMeshGPUDeformer.h"

#define GLH_EXT_SINGLE_FILE
#include "glh_extensions.h"

#ifndef WIN32
  #define GL_GLEXT_PROTOTYPES 1
#endif

#if defined(linux)
  GLAPI void APIENTRY glDeleteBuffersARB (GLsizei, const GLuint *);
  GLAPI void APIENTRY glBindBufferARB (GLenum, GLuint);
  GLAPI void APIENTRY glGenBuffersARB (GLsizei, GLuint *);
  GLAPI void APIENTRY glBufferDataARB (GLenum, GLsizeiptrARB, const GLvoid *, GLenum);
#endif

#ifdef OBJMESHGPUDEFORMER_USING_VBOS
  #include "vbo.h"
#endif

ObjMeshGPUDeformer::ObjMeshGPUDeformer()
{
  numGroupTriangles = NULL;

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    vboID = NULL;
    vboEdgesID = NULL;
  #endif

  vboNormalEnabled = NULL;
  vboTex1Enabled = NULL;
}

void ObjMeshGPUDeformer::Init(ObjMesh * mesh_, ObjMeshRender * meshRender_, int renderingMode_) 
{
  mesh = mesh_;
  meshRender = meshRender_;
  renderingMode = renderingMode_;

  printf("Initializing ObjMeshGPUDeformer...\n");

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    // init the GL_ARB_vertex_buffer_object extension
    if (!InitializeVBOs())
    {
      printf("Unable to load the GL_ARB_vertex_buffer_object extension\n");
      throw 31;
    }
    printf("Detected extension: GL_ARB_vertex_buffer_object\n");
  #endif

  // default lighting: single light of intensity 1 at the origin
  // set the other lights away from the object 
  // (to avoid a vertex accidentally coinciding with a light source (which will cause a division by zero in the shader when normalizing the light vector))
  memset(lightPos, 0, sizeof(float) * 16);
  memset(lightIntensity, 0, sizeof(float) * 4);
  lightIntensity[0] = 1.0;
  ambientIntensity = 0.0;

  // set the other three lights to be far away and of zero intensity
  float vtx_max[4] = { FLT_MAX, FLT_MAX, FLT_MAX, 1.0 };

  SetLightPosition(1, vtx_max);
  SetLightPosition(2, vtx_max);
  SetLightPosition(3, vtx_max);

  ambientEnabled = true;
  diffuseEnabled = true;
  specularEnabled = true;

  numVertices = mesh->getNumVertices();
  numGroups = mesh->getNumGroups();
  numGroupTriangles = (int*) malloc (sizeof(int) * numGroups);

  // count num mesh triangles
  numTriangles = 0;
  for(int groupNo=0; groupNo<numGroups; groupNo++)
  {
    const ObjMesh::Group * groupHandle = mesh->getGroupHandle(groupNo);

    // count num triangles
    numGroupTriangles[groupNo] = 0;
    for(unsigned int iFace = 0; iFace < groupHandle->getNumFaces(); iFace++)
    {
      const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(iFace);
      numGroupTriangles[groupNo] += faceHandle->getNumVertices() - 2;
    }
    numTriangles += numGroupTriangles[groupNo];
  }

  InitVertexDeformationTexture();

  printf ("ObjMeshGPUDeformer initialization complete.\n");
}

ObjMeshGPUDeformer::~ObjMeshGPUDeformer()
{
  free(numGroupTriangles);

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    if (vboID != NULL)
    {
      glDeleteBuffersARB(5 * numGroups, vboID);
      free(vboID);
    }
    if (vboEdgesID != NULL)
    {
      glDeleteBuffersARB(3 * numGroups, vboEdgesID);
      free(vboEdgesID);
    }
  #endif

  free(vboNormalEnabled);
  free(vboTex1Enabled);
}

void ObjMeshGPUDeformer::InitVertexDeformationTexture()
{
  gpgpuVertexTextureCoordinates = (float*) calloc (2*numVertices, sizeof(float));
  vertexDeformationTextureSize = 1 << (1+(int)(0.5*log(1.0f * numVertices)/log (2.0f)));

  int row = 0;
  int column = 0;
  for(int i=0; i< numVertices; i++)
  {
    if (column >= vertexDeformationTextureSize)
    {
      column = 0;
      row++;
    }

    gpgpuVertexTextureCoordinates[2*i+0] = ((0.5 + column) / vertexDeformationTextureSize);
    gpgpuVertexTextureCoordinates[2*i+1] = ((0.5 + row) / vertexDeformationTextureSize);

    column++;
  }

  printf("Vertex texture coordinates initialized. Size of vertex texture is %d x %d .\n", vertexDeformationTextureSize, vertexDeformationTextureSize);

  glEnable(GL_TEXTURE_2D);

  // Create the dynamic texture (render-to-texture target)

  // to map it to 2nd texture unit:
  //glActiveTextureARB(GL_TEXTURE1_ARB); 

  GLenum err;
  if ((err = glGetError()) != GL_NO_ERROR)
  {
    const GLubyte * errString = gluErrorString(err);
    printf("Warning at (1): %s \n",errString);
  }

  glGenTextures( 1, &vertexDeformationTextureID );
  glBindTexture( GL_TEXTURE_2D, vertexDeformationTextureID );

  glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST );

/*
  float * test = (float*) malloc (sizeof(float) * vertexDeformationTextureSize * vertexDeformationTextureSize * 4);
  for(int i1 = 0 ; i1 < vertexDeformationTextureSize; i1++)
    for(int j1 = 0 ; j1 < vertexDeformationTextureSize; j1++)
    {
      test[(j1 * vertexDeformationTextureSize + i1) * 4 + 0] = 0.0;
      test[(j1 * vertexDeformationTextureSize + i1) * 4 + 1] = 1.0;
      test[(j1 * vertexDeformationTextureSize + i1) * 4 + 2] = 0.0;
      test[(j1 * vertexDeformationTextureSize + i1) * 4 + 3] = 1.0;
    }
  free(test);

  if ((err = glGetError()) != GL_NO_ERROR)
  {
    const GLubyte * errString = gluErrorString(err);
    printf("Warning at (2): %s \n",errString);
  }
*/

  float * zeroTexture = (float*) calloc (vertexDeformationTextureSize * vertexDeformationTextureSize * 4, sizeof(float));

  //int internalFormat = GL_RGBA_FLOAT16_ATI;
  //int internalFormat = GL_RGBA_FLOAT32_ATI;
  //int internalFormat = GL_RGBA16F_ARB;
  int internalFormat = GL_RGBA32F_ARB;
  glTexImage2D( GL_TEXTURE_2D, 0, internalFormat,
    vertexDeformationTextureSize, vertexDeformationTextureSize, 0, 
    GL_RGBA, GL_FLOAT, zeroTexture );

  free(zeroTexture);

  if ((err = glGetError()) != GL_NO_ERROR)
  {
    const GLubyte * errString = gluErrorString(err);
    printf("Warning at (3): %s \n", errString);
  }
}

int ObjMeshGPUDeformer::glh_extension_supported(const char *extension)
{
    static const GLubyte *extensions = NULL;
    const GLubyte *start;
    GLubyte *where, *terminator;
    
    // Extension names should not have spaces.
    where = (GLubyte *) strchr(extension, ' ');
    if (where || *extension == '\0')
      return 0;
    
    if (!extensions)
      extensions = glGetString(GL_EXTENSIONS);

    //printf("%s\n", extensions);

    // It takes a bit of care to be fool-proof about parsing the
    // OpenGL extensions string.  Don't be fooled by sub-strings,
    // etc.
    start = extensions;
    for (;;) 
    {
        where = (GLubyte *) strstr((const char *) start, extension);
        if (!where)
            break;
        terminator = where + strlen(extension);
        if (where == start || *(where - 1) == ' ') 
        {
            if (*terminator == ' ' || *terminator == '\0') 
            {
                return 1;
            }
        }
        start = terminator;
    }
    return 0;
}

void ObjMeshGPUDeformer::SetLightPosition(int lightID, float pos[4])
{
  memcpy(&(lightPos[4*lightID]), pos, sizeof(float) * 4);
}

void ObjMeshGPUDeformer::SetLightIntensity(int lightID, float intensity)
{ 
  lightIntensity[lightID] = intensity;
}

void ObjMeshGPUDeformer::SetAmbientIntensity(float intensity)
{ 
  ambientIntensity = intensity;
}

void ObjMeshGPUDeformer::RenderMaster(int masterMode, void * data)
{
  glActiveTextureARB(GL_TEXTURE0_ARB); 
  glDisable(GL_TEXTURE_2D);

  cgGLEnableProfile(VertexProfile);

  if (masterMode == 0)
  {
    cgGLBindProgram(VertexPass2Program);

    if (renderingMode & OBJMESHRENDER_TEXTURE)
    {
      cgGLEnableProfile(FragmentProfile);
      cgGLBindProgram(FragmentPass2Program);
    }

    cgGLSetStateMatrixParameter(ModelViewProjParam,
      CG_GL_MODELVIEW_PROJECTION_MATRIX,
      CG_GL_MATRIX_IDENTITY);

    cgGLSetStateMatrixParameter(ModelViewITParam,
      CG_GL_MODELVIEW_MATRIX,
      CG_GL_MATRIX_INVERSE_TRANSPOSE);

    cgGLSetParameter4f(LightPos1Param, lightPos[0], lightPos[1], lightPos[2], 1);  
    cgGLSetParameter4f(LightPos2Param, lightPos[4], lightPos[5], lightPos[6], 1);  
    cgGLSetParameter4f(LightPos3Param, lightPos[8], lightPos[9], lightPos[10], 1);  
    cgGLSetParameter4f(LightPos4Param, lightPos[12], lightPos[13], lightPos[14], 1); 
    cgGLSetParameter1f(Light1IntensityParam, lightIntensity[0]);  
    cgGLSetParameter1f(Light2IntensityParam, lightIntensity[1]);  
    cgGLSetParameter1f(Light3IntensityParam, lightIntensity[2]);  
    cgGLSetParameter1f(Light4IntensityParam, lightIntensity[3]);  

    cgGLSetParameter1f(AmbientIntensityParam, ambientIntensity);  
  }

  if (masterMode == 1)
  {
    cgGLBindProgram(VertexPass2ProgramShadow);

    cgGLSetStateMatrixParameter(ModelViewProjShadowParam,
      CG_GL_MODELVIEW_PROJECTION_MATRIX,
      CG_GL_MATRIX_IDENTITY);

    double intensity = *(double*)data;
    cgGLSetParameter1f(ShadowIntensityParam, intensity);
  }

  if (renderingMode & OBJMESHRENDER_TEXTURE)
  {
    cgGLEnableProfile(FragmentProfile);
    cgGLBindProgram(FragmentPass2Program);
  }

  cgGLSetStateMatrixParameter(ModelViewProjParam,
    CG_GL_MODELVIEW_PROJECTION_MATRIX,
    CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(ModelViewITParam,
    CG_GL_MODELVIEW_MATRIX,
    CG_GL_MATRIX_INVERSE_TRANSPOSE);

  cgGLSetParameter4f(LightPos1Param, lightPos[0], lightPos[1], lightPos[2], 1);  
  cgGLSetParameter4f(LightPos2Param, lightPos[4], lightPos[5], lightPos[6], 1);  
  cgGLSetParameter4f(LightPos3Param, lightPos[8], lightPos[9], lightPos[10], 1);  
  cgGLSetParameter4f(LightPos4Param, lightPos[12], lightPos[13], lightPos[14], 1); 
  cgGLSetParameter1f(Light1IntensityParam, lightIntensity[0]);  
  cgGLSetParameter1f(Light2IntensityParam, lightIntensity[1]);  
  cgGLSetParameter1f(Light3IntensityParam, lightIntensity[2]);  
  cgGLSetParameter1f(Light4IntensityParam, lightIntensity[3]);  

  cgGLSetParameter1f(AmbientIntensityParam, ambientIntensity);  

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);
  BindRT();

  for(int groupNo=0; groupNo<numGroups; groupNo++)
  {
    const ObjMesh::Group * groupHandle = mesh->getGroupHandle(groupNo);
    int materialIndex = groupHandle->getMaterialIndex();
    const ObjMesh::Material * materialHandle = mesh->getMaterialHandle(materialIndex);

    float alpha = (float)(materialHandle->getAlpha());

    Vec3d Ka = materialHandle->getKa();
    float Kav[4] = { (float)Ka[0], (float)Ka[1], (float)Ka[2], alpha };

    Vec3d Kd = materialHandle->getKd();
    float Kdv[4] = { (float)Kd[0], (float)Kd[1], (float)Kd[2], alpha };

    Vec3d Ks = materialHandle->getKs();
    float Ksv[4] = { (float)Ks[0], (float)Ks[1], (float)Ks[2], alpha };

    float shininess = materialHandle->getShininess();

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Kav);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kdv);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ksv);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    float zero[4] = { 0, 0, 0, 0 };
    cgGLSetParameter4fv(KaParam, ambientEnabled ? Kav : zero);
    cgGLSetParameter4fv(KdParam, diffuseEnabled ? Kdv : zero);
    cgGLSetParameter4fv(KsParam, specularEnabled ? Ksv : zero);
    cgGLSetParameter1f(shininessParam, shininess);

    if ((renderingMode & OBJMESHRENDER_TEXTURE) && (materialHandle->hasTextureFilename()))
    {
      ObjMeshRender::Texture * textureHandle = meshRender->getTextureHandle(materialIndex);
      GLuint tname = textureHandle->getTexture();
      glActiveTextureARB(GL_TEXTURE1_ARB); 
      glBindTexture(GL_TEXTURE_2D, tname);
      glEnable(GL_TEXTURE_2D);
    }

    #ifdef OBJMESHGPUDEFORMER_USING_VBOS
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+0]);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, 0, 0);

      if (vboNormalEnabled[groupNo])
      {
        glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+1]);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, 0);
      }

      glClientActiveTextureARB(GL_TEXTURE0_ARB);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+2]);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      glTexCoordPointer(2, GL_FLOAT, 0, 0);

      if (vboTex1Enabled[groupNo])
      {
        glClientActiveTextureARB(GL_TEXTURE1_ARB);
        glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+3]);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, 0, 0);
      }

      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vboID[5*groupNo+4]);
      glDrawElements(GL_TRIANGLES, 3 * numGroupTriangles[groupNo], GL_UNSIGNED_INT, 0);

      // turn off VBOs
      glDisableClientState(GL_VERTEX_ARRAY);
      if (vboNormalEnabled[groupNo])
        glDisableClientState(GL_NORMAL_ARRAY);
      glClientActiveTextureARB(GL_TEXTURE0_ARB);
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);
      if (vboTex1Enabled[groupNo])
      {
        glClientActiveTextureARB(GL_TEXTURE1_ARB);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
      }
      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

    #else
      glCallList(displayListStart + groupNo);
    #endif

    if ((masterMode == 0) && (renderingMode & OBJMESHRENDER_TEXTURE) && (materialHandle->hasTextureFilename()))
    {
      glDisable(GL_TEXTURE_2D);
      glActiveTextureARB(GL_TEXTURE0_ARB); 
    }
  }

  UnbindRT();

  if ((masterMode == 0) && (renderingMode & OBJMESHRENDER_TEXTURE))
    cgGLDisableProfile(FragmentProfile);

  cgGLDisableProfile(VertexProfile);
}

void ObjMeshGPUDeformer::Render()
{
  RenderMaster(0);
}

void ObjMeshGPUDeformer::RenderShadow(double intensity)
{
  RenderMaster(1, &intensity);
}

void ObjMeshGPUDeformer::RenderVertices()
{
  glDisable(GL_TEXTURE_2D);

  cgGLEnableProfile(VertexProfile);
  cgGLBindProgram(VertexPass2ProgramPoints);

  cgGLSetStateMatrixParameter(ModelViewProjPointParam,
    CG_GL_MODELVIEW_PROJECTION_MATRIX,
    CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(ModelViewITPointParam,
    CG_GL_MODELVIEW_MATRIX,
    CG_GL_MATRIX_INVERSE_TRANSPOSE);

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);
  BindRT();
    glCallList(displayListPoints);
  UnbindRT();

  cgGLDisableProfile(VertexProfile);
}

void ObjMeshGPUDeformer::RenderEdges()
{
  glDisable(GL_TEXTURE_2D);

  cgGLEnableProfile(VertexProfile);
  cgGLBindProgram(VertexPass2ProgramEdges);

  cgGLSetStateMatrixParameter(ModelViewProjEdgeParam,
    CG_GL_MODELVIEW_PROJECTION_MATRIX,
    CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(ModelViewITEdgeParam,
    CG_GL_MODELVIEW_MATRIX,
    CG_GL_MATRIX_INVERSE_TRANSPOSE);

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);
  BindRT();

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    glEnableClientState(GL_VERTEX_ARRAY);
    glClientActiveTextureARB(GL_TEXTURE0_ARB);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    for(int groupNo=0; groupNo<numGroups; groupNo++)
    {
      // render via VBOs
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+0]);
      glVertexPointer(3, GL_FLOAT, 0, 0);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+1]);
      glTexCoordPointer(2, GL_FLOAT, 0, 0);

      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+2]);

      int numGroupEdges = 3 * numGroupTriangles[groupNo];
      glDrawElements(GL_LINES, 2 * numGroupEdges, GL_UNSIGNED_INT, 0);
    }

    // unbind VBOs
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
  #else
    for(int groupNo=0; groupNo<numGroups; groupNo++)
      glCallList(displayListEdgesStart + groupNo);
  #endif

  UnbindRT();

  cgGLDisableProfile(VertexProfile);
}

void ObjMeshGPUDeformer::MakeDisplayLists(int mode)
{
  PrintGLerror("before making display lists");
  MakeDisplayListsTriangles(mode);
  MakeDisplayListsPoints();
  MakeDisplayListsEdges();
  free(gpgpuVertexTextureCoordinates);
  gpgpuVertexTextureCoordinates = NULL;
}

void ObjMeshGPUDeformer::MakeDisplayListsTriangles(int mode)
{
  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    printf("Creating VBOs...\n");
  #else
    printf("Creating display lists...\n");
  #endif

  // mode must be either OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL
  //                  or OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL | OBJMESHRENDER_TEXTURE

  if ((mode != (OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL)) &&
    (mode != (OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL | OBJMESHRENDER_TEXTURE)))
  {
    printf("Error: invalid rendering mode.\n");
    exit(1);
  }

  if (mode & OBJMESHRENDER_COLOR)
    glEnable(GL_COLOR_MATERIAL);
  else if (mode & OBJMESHRENDER_MATERIAL)
    glDisable(GL_COLOR_MATERIAL);

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    vboID = (GLuint*) malloc (sizeof(GLuint) * 5 * numGroups);
    glGenBuffersARB(5 * numGroups, vboID);
  #else
    // create standard display lists
    displayListStart = glGenLists(numGroups);
  #endif

  // create the array data
  vboNormalEnabled = (int*) calloc (numGroups, sizeof(int));
  vboTex1Enabled = (int*) calloc (numGroups, sizeof(int));

  for(int groupNo=0; groupNo<numGroups; groupNo++)
  {
    // allocate space for group vbo data
    float * vtxBuffer = (float*) malloc (sizeof(float) * 9 * numGroupTriangles[groupNo]);
    float * normalBuffer = NULL;
    if (mode & OBJMESHRENDER_SMOOTH)
    {
      vboNormalEnabled[groupNo] = 1;
      normalBuffer = (float*) malloc (sizeof(float) * 9 * numGroupTriangles[groupNo]);
    }
    float * tex0Buffer = (float*) malloc (sizeof(float) * 6 * numGroupTriangles[groupNo]);
    float * tex1Buffer = NULL;
    if (mode & OBJMESHRENDER_TEXTURE)
    {
      vboTex1Enabled[groupNo] = 1;
      tex1Buffer = (float*) malloc (sizeof(float) * 6 * numGroupTriangles[groupNo]);
    }
    GLuint * indexBuffer = (GLuint*) malloc (sizeof(GLuint) * 3 * numGroupTriangles[groupNo]);

    const ObjMesh::Group * groupHandle = mesh->getGroupHandle(groupNo);

    int triangleCount = 0;
    for(unsigned int iFace = 0; iFace < groupHandle->getNumFaces(); iFace++)
    {
      const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(iFace);

      // triangulate the face on the fly
      for(unsigned int iVtx = 0; iVtx < faceHandle->getNumVertices() - 2; iVtx++)
      {
        unsigned int triangleVertex[3] = { 0, iVtx + 1, iVtx + 2 };

        for (int vtx=0; vtx<3; vtx++)
        {
          const ObjMesh::Vertex * vertex = faceHandle->getVertexHandle(triangleVertex[vtx]);
          Vec3d pos = mesh->getPosition(*vertex);
          for (int dof=0; dof<3; dof++)
            vtxBuffer[9*triangleCount + 3*vtx + dof] = pos[dof];

          if (vboNormalEnabled[groupNo])
          {
            Vec3d normal = mesh->getNormal(*vertex);
            for (int dof=0; dof<3; dof++)
              normalBuffer[9*triangleCount + 3*vtx + dof] = normal[dof];
          }

          int vertexIndex = vertex->getPositionIndex();
          Vec3d uv = Vec3d(gpgpuVertexTextureCoordinates[2*vertexIndex+0], 
                           gpgpuVertexTextureCoordinates[2*vertexIndex+1], 0);
          for (int dof=0; dof<2; dof++)
            tex0Buffer[6*triangleCount + 2*vtx + dof] = uv[dof];

          if (vboTex1Enabled[groupNo])
          {
            Vec3d uv = mesh->getTextureCoordinate(*vertex);
            for (int dof=0; dof<2; dof++)
              tex1Buffer[6*triangleCount + 2*vtx + dof] = uv[dof];
          }
        }

        for (int dof=0; dof<3; dof++)
          indexBuffer[3*triangleCount+dof] = 3*triangleCount+dof;

        triangleCount++;
      }
    }

    #ifdef OBJMESHGPUDEFORMER_USING_VBOS
      // upload the VBOs
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+0]);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 9 * numGroupTriangles[groupNo], vtxBuffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+1]);
      if (vboNormalEnabled[groupNo])
        glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 9 * numGroupTriangles[groupNo], normalBuffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+2]);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 6 * numGroupTriangles[groupNo], tex0Buffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboID[5*groupNo+3]);
      if (vboTex1Enabled[groupNo])
        glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 6 * numGroupTriangles[groupNo], tex1Buffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vboID[5*groupNo+4]); 
      glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLuint) * 3 * numGroupTriangles[groupNo], indexBuffer, GL_STATIC_DRAW_ARB);

    #else

      PrintGLerror("before starting a new list");
      glNewList(displayListStart + groupNo, GL_COMPILE);

      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, 0, vtxBuffer);

      if (vboNormalEnabled[groupNo])
      {
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, normalBuffer);
      }

      glClientActiveTextureARB(GL_TEXTURE0_ARB); 
      glTexCoordPointer(2, GL_FLOAT, 0, tex0Buffer); 
      glEnableClientState(GL_TEXTURE_COORD_ARRAY); 

      if (vboTex1Enabled[groupNo])
      {
        glClientActiveTextureARB(GL_TEXTURE1_ARB); 
        glTexCoordPointer(2, GL_FLOAT, 0, tex1Buffer); 
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      }

      glDrawElements(GL_TRIANGLES, 3 * numGroupTriangles[groupNo], GL_INT, indexBuffer);

      glDisableClientState(GL_VERTEX_ARRAY);
      glDisableClientState(GL_NORMAL_ARRAY);
      glClientActiveTextureARB(GL_TEXTURE0_ARB); 
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);
      glClientActiveTextureARB(GL_TEXTURE1_ARB); 
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);

      glEndList();    

      char msg[4096];
      sprintf(msg,"displayList creation, groupNo: %d\n",groupNo);
      PrintGLerror(msg);
    #endif

    // de-allocate space for group vbo data
    free(vtxBuffer);
    free(normalBuffer);
    free(tex0Buffer);
    free(tex1Buffer);
    free(indexBuffer);
  }

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    // unbind buffers
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  #endif
}

void ObjMeshGPUDeformer::MakeDisplayListsPoints()
{
  printf("Creating display list for points...\n");

  displayListPoints = glGenLists(1);

  glNewList(displayListPoints, GL_COMPILE);
  glBegin(GL_POINTS);
  for (int i=0; i < (int)mesh->getNumVertices(); i++)
  {
    float s = gpgpuVertexTextureCoordinates[2*i+0];
    float t = gpgpuVertexTextureCoordinates[2*i+1];
    glMultiTexCoord2fARB(GL_TEXTURE0_ARB, s, t);

    Vec3d pos = mesh->getPosition(i);
    glVertex3f(pos[0], pos[1], pos[2]);
  }
  glEnd();
  glEndList();
}    

void ObjMeshGPUDeformer::MakeDisplayListsEdges()
{
  printf("Creating display list for edges...\n");

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    vboEdgesID = (GLuint*) malloc (sizeof(GLuint) * 3 * numGroups);
    glGenBuffersARB(3 * numGroups, vboEdgesID);
  #else
    displayListEdgesStart = glGenLists(numGroups);
  #endif

  for(int groupNo=0; groupNo<numGroups; groupNo++)
  {
    int numGroupEdges = 3 * numGroupTriangles[groupNo];
    float * vtxBuffer = (float*) malloc (sizeof(float) * 6 * numGroupEdges);
    float * stBuffer = (float*) malloc (sizeof(float) * 4 * numGroupEdges);
    GLuint * indexBuffer = (GLuint*) malloc (sizeof(GLuint) * 2 * numGroupEdges);

    const ObjMesh::Group * groupHandle = mesh->getGroupHandle(groupNo);
    int edgeCount = 0;

    for(unsigned int iFace = 0; iFace < groupHandle->getNumFaces(); iFace++)
    {
      const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(iFace);

      for(unsigned int iVtx = 0; iVtx < faceHandle->getNumVertices() - 1; iVtx++)
      {
        unsigned int edgeVertex[2] = { iVtx, iVtx + 1 };

        for (int vtx=0; vtx<2; vtx++)
        {
          const ObjMesh::Vertex * vertex = faceHandle->getVertexHandle(edgeVertex[vtx]);
          int vertexIndex = vertex->getPositionIndex();
          Vec3d pos = mesh->getPosition(*vertex);
          for (int dof=0; dof<3; dof++)
            vtxBuffer[6*edgeCount + 3*vtx + dof] = pos[dof];

          float s = gpgpuVertexTextureCoordinates[2*vertexIndex+0];
          float t = gpgpuVertexTextureCoordinates[2*vertexIndex+1];
          float st[2] = {s, t};

          for (int dof=0; dof<2; dof++)
            stBuffer[4*edgeCount + 2*vtx + dof] = st[dof];

          for (int dof=0; dof<2; dof++)
            indexBuffer[2 * edgeCount + dof] = 2 * edgeCount + dof;
        }
       
        edgeCount++;
      }
    }

    #ifdef OBJMESHGPUDEFORMER_USING_VBOS
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+0]);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 6 * numGroupEdges, vtxBuffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+1]);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 4 * numGroupEdges, stBuffer, GL_STATIC_DRAW_ARB);

      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vboEdgesID[3*groupNo+2]); 
      glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLuint) * 2 * numGroupEdges, indexBuffer, GL_STATIC_DRAW_ARB);
    #else
      PrintGLerror("before starting a new edge list");
      glNewList(displayListEdgesStart + groupNo, GL_COMPILE);

      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, 0, vtxBuffer);

      glClientActiveTextureARB(GL_TEXTURE0_ARB); 
      glTexCoordPointer(2, GL_FLOAT, 0, stBuffer); 
      glEnableClientState(GL_TEXTURE_COORD_ARRAY); 

      glDrawElements(GL_LINES, 2 * numGroupEdges, GL_INT, indexBuffer);

      glDisableClientState(GL_VERTEX_ARRAY);
      glClientActiveTextureARB(GL_TEXTURE0_ARB); 
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);

      glEndList();    
    #endif

    free(indexBuffer);
    free(stBuffer);
    free(vtxBuffer);
  }

  #ifdef OBJMESHGPUDEFORMER_USING_VBOS
    // unbind buffer
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  #endif
}    

void ObjMeshGPUDeformer::ReadBack_u(float * u)
{
  PrintGLerror("Entering readback");
  float * buffer = (float*) malloc (sizeof(float) * 4 * vertexDeformationTextureSize * vertexDeformationTextureSize);

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);
  glGetTexImage(GL_TEXTURE_2D,0,GL_RGBA,GL_FLOAT,buffer);

  for(int i=0; i< numVertices; i++)
  {
    u[3*i+0] = buffer[4*i+0];
    u[3*i+1] = buffer[4*i+1];
    u[3*i+2] = buffer[4*i+2];
  }
  free(buffer);
  PrintGLerror("Leaving readback");
}

void ObjMeshGPUDeformer::ReadBack_u(double * u)
{
  int numVertices3 = 3 * numVertices;
  float * buffer = (float*) malloc (sizeof(float) * numVertices3);
  ReadBack_u(buffer);
  for(int i=0; i<numVertices3; i++)
    u[i] = buffer[i];
  free(buffer);
}

void ObjMeshGPUDeformer::DeleteCGShaders()
{
  if (VertexPass2Program) 
    cgDestroyProgram(VertexPass2Program);
  
  if (VertexPass2ProgramShadow)
    cgDestroyProgram(VertexPass2ProgramShadow);

  //if (VertexPass2ProgramDeformedNormals)
    //cgDestroyProgram(VertexPass2ProgramDeformedNormals);

  if (VertexPass2ProgramPoints) 
    cgDestroyProgram(VertexPass2ProgramPoints);
    
  if (VertexPass2ProgramEdges)
    cgDestroyProgram(VertexPass2ProgramEdges);

  if (FragmentPass2Program)
    cgDestroyProgram(FragmentPass2Program);
    
  if (Context)
    cgDestroyContext(Context);
}

void ObjMeshGPUDeformer::cgErrorCallback(void)
{
  CGerror LastError = cgGetError();

  if (LastError){
    switch(LastError){
      case CG_COMPILER_ERROR:printf("compiler error\n");break;
      case CG_INVALID_PARAMETER_ERROR:printf("invalid parameter error\n");break;
      case CG_INVALID_PROFILE_ERROR:printf("invalid profile error\n");break;
      case CG_INVALID_VALUE_TYPE_ERROR:printf("invalid value type error\n");break;
      case CG_NOT_MATRIX_PARAM_ERROR:printf("not matrix param error\n");break;
      case CG_INVALID_ENUMERANT_ERROR:printf("invalid enumerant error\n");break;
      case CG_NOT_4x4_MATRIX_ERROR:printf("not 4x4 matrix error\n");break;
      case CG_FILE_READ_ERROR:printf("file read error\n");break;
      case CG_FILE_WRITE_ERROR:printf("file write error\n");break;
      case CG_MEMORY_ALLOC_ERROR:printf("memory alloc error\n");break;
      case CG_INVALID_CONTEXT_HANDLE_ERROR:printf("invalid context handle error\n");break;
      case CG_INVALID_PROGRAM_HANDLE_ERROR:printf("invalid program handle error\n");break;
      case CG_INVALID_PARAM_HANDLE_ERROR:printf("invalid param handle error\n");break;
      case CG_UNKNOWN_PROFILE_ERROR:printf("unknown profile error\n");break;
      case CG_VAR_ARG_ERROR:printf("var arg error\n");break;
      case CG_INVALID_DIMENSION_ERROR:printf("invalid dimension error\n");break;
      case CG_ARRAY_PARAM_ERROR:printf("array param error\n");break;
      case CG_OUT_OF_ARRAY_BOUNDS_ERROR:printf("out of array bounds error\n");break;
      case CG_NO_ERROR:printf("no error\n");break;
      case CG_PROGRAM_LOAD_ERROR:printf("program load error\n");break;
      case CG_PROGRAM_BIND_ERROR:printf("program bind error\n");break;
      case CG_PROGRAM_NOT_LOADED_ERROR:printf("program not loaded error\n");break;
      case CG_UNSUPPORTED_GL_EXTENSION_ERROR:printf("unsupported gl extension error\n");break;
      case CG_NVPARSE_ERROR:printf("nvparse error\n");break;
      default:printf("unknown error\n");
    }

  }

  if(LastError)
  {
    const char *Listing = cgGetLastListing(Context);
    printf("\n---------------------------------------------------\n");
    printf("%s\n\n", cgGetErrorString(LastError));
    printf("%s\n", Listing);
    printf("---------------------------------------------------\n");
    printf("Cg error, exiting...\n");
    exit(0);
  }
}

void ObjMeshGPUDeformer::PrintGLerror( const char *msg )
{
  GLenum errCode;
  const GLubyte *errStr;

  if ((errCode = glGetError()) != GL_NO_ERROR) 
  {
    errStr = gluErrorString(errCode);
    printf("OpenGL ERROR: %s: %s\n", errStr, msg);
    exit(1);
  }
}

#ifdef WIN32
  #include <malloc.h>

  void ObjMeshGPUDeformer::heap_check_()
  {
    int heapstatus = _heapchk();
    switch( heapstatus )
    {
      case _HEAPOK:
        printf(" OK - heap is fine\n" );
      break;
      case _HEAPEMPTY:
        printf(" OK - heap is empty\n" );
      break;
      case _HEAPBADBEGIN:
        printf( "ERROR - bad start of heap\n" );
      break;
      case _HEAPBADNODE:
        printf( "ERROR - bad node in heap\n" );
      break;
    }
  }
#endif

// static members
CGcontext ObjMeshGPUDeformer::Context; 
CGprofile ObjMeshGPUDeformer::VertexProfile;
CGprogram ObjMeshGPUDeformer::VertexPass2Program;
CGprogram ObjMeshGPUDeformer::VertexPass2ProgramShadow;
CGprogram ObjMeshGPUDeformer::VertexPass2ProgramPoints;
CGprogram ObjMeshGPUDeformer::VertexPass2ProgramEdges;
//CGprogram ObjMeshGPUDeformer::VertexPass2ProgramDeformedNormals;
CGprofile ObjMeshGPUDeformer::FragmentProfile;
CGprogram ObjMeshGPUDeformer::FragmentPass2Program;

// vertex shader parameters
CGparameter ObjMeshGPUDeformer::texParam;
CGparameter ObjMeshGPUDeformer::fragmentTexParam;
CGparameter ObjMeshGPUDeformer::KaParam;
CGparameter ObjMeshGPUDeformer::KdParam;
CGparameter ObjMeshGPUDeformer::KsParam;
CGparameter ObjMeshGPUDeformer::shininessParam;
CGparameter ObjMeshGPUDeformer::ModelViewProjParam;
CGparameter ObjMeshGPUDeformer::ModelViewITParam;
CGparameter ObjMeshGPUDeformer::LightPos1Param;
CGparameter ObjMeshGPUDeformer::LightPos2Param;
CGparameter ObjMeshGPUDeformer::LightPos3Param;
CGparameter ObjMeshGPUDeformer::LightPos4Param;
CGparameter ObjMeshGPUDeformer::Light1IntensityParam;
CGparameter ObjMeshGPUDeformer::Light2IntensityParam;
CGparameter ObjMeshGPUDeformer::Light3IntensityParam;
CGparameter ObjMeshGPUDeformer::Light4IntensityParam;
CGparameter ObjMeshGPUDeformer::AmbientIntensityParam;
CGparameter ObjMeshGPUDeformer::vertexTextureCoordParam;

// vertex shader shadow parameters
CGparameter ObjMeshGPUDeformer::ModelViewProjShadowParam;
CGparameter ObjMeshGPUDeformer::ShadowIntensityParam;

// vertex shader for points parameters
CGparameter ObjMeshGPUDeformer::ModelViewProjPointParam;
CGparameter ObjMeshGPUDeformer::ModelViewITPointParam;

// vertex shader for edges parameters
CGparameter ObjMeshGPUDeformer::ModelViewProjEdgeParam;
CGparameter ObjMeshGPUDeformer::ModelViewITEdgeParam;
 
// pass 2 fragment shader parameters
CGparameter ObjMeshGPUDeformer::texPass2Param;

/*
void ObjMeshGPUDeformer::RenderWithDeformedNormals()
{
  glActiveTextureARB(GL_TEXTURE0_ARB); 
  glDisable(GL_TEXTURE_2D);

  cgGLEnableProfile(VertexProfile);

  cgGLBindProgram(VertexProgramDeformedNormals);

  if (renderingMode & OBJMESHRENDER_TEXTURE)
  {
    cgGLEnableProfile(FragmentProfile);
    cgGLBindProgram(FragmentPass2Program);
  }

  cgGLSetStateMatrixParameter(ModelViewProjDeformedNormalsParam,
    CG_GL_MODELVIEW_PROJECTION_MATRIX,
    CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(ModelViewITDeformedNormalsParam,
    CG_GL_MODELVIEW_MATRIX,
    CG_GL_MATRIX_INVERSE_TRANSPOSE);

  cgGLSetParameter4f(LightPos1DeformedNormalsParam, lightPos[0], lightPos[1], lightPos[2], 1);  
  cgGLSetParameter4f(LightPos2DeformedNormalsParam, lightPos[4], lightPos[5], lightPos[6], 1);  
  cgGLSetParameter4f(LightPos3DeformedNormalsParam, lightPos[8], lightPos[9], lightPos[10], 1);  
  cgGLSetParameter4f(LightPos4DeformedNormalsParam, lightPos[12], lightPos[13], lightPos[14], 1); 

  cgGLSetParameter1f(Light1IntensityDeformedNormalsParam, lightIntensity[0]);  
  cgGLSetParameter1f(Light2IntensityDeformedNormalsParam, lightIntensity[1]);  
  cgGLSetParameter1f(Light3IntensityDeformedNormalsParam, lightIntensity[2]);  
  cgGLSetParameter1f(Light4IntensityDeformedNormalsParam, lightIntensity[3]);  

  cgGLSetParameter1f(AmbientIntensityDeformedNormalsParam, ambientIntensity);  

  glBindTexture(GL_TEXTURE_2D, vertexDeformationTextureID);

  if( !wglBindTexImageARB( pbuffer.hPBuffer, WGL_FRONT_LEFT_ARB ) )
  {
    printf("Could not bind p-buffer to render texture!");
    exit(-1);
  }

  GLMgroup* group;
  GLMmaterial* material;

  int i=0;
  group = model->groups;
  while (group) 
  {
    material = &model->materials[group->material];
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material->ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material->diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material->specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material->shininess);

    cgGLSetParameter4fv(KaParam, material->ambient);
    cgGLSetParameter4fv(KdParam, material->diffuse);
    //printf("%d | %G %G %G %G | %G\n", group->numtriangles,
      //material->specular[0], material->specular[1], 
      //material->specular[2], material->specular[3], material->shininess);
    cgGLSetParameter4fv(KsParam, material->specular);
    cgGLSetParameter1f(shininessParam, material->shininess);

    if ((renderingMode & OBJMESHRENDER_TEXTURE) && (material->textureData != NULL))
    {
      GLuint tname = material->textureName;
      glActiveTextureARB(GL_TEXTURE1_ARB); 
      glBindTexture(GL_TEXTURE_2D, tname);
      glEnable(GL_TEXTURE_2D);
    }

    glCallList(displayListStart+i);

    if ((renderingMode & OBJMESHRENDER_TEXTURE) && (material->textureData != NULL))
    {
      glDisable(GL_TEXTURE_2D);
      glActiveTextureARB(GL_TEXTURE0_ARB); 
    }

    group = group->next;
    i++;
  }

  // Before we forget, we need to make sure that the p-buffer has been 
  // released from the dynamic "render-to" texture.
  if( !wglReleaseTexImageARB( pbuffer.hPBuffer, WGL_FRONT_LEFT_ARB ) )
  {
    printf("Could not release p-buffer from render texture!");
    exit(-1);
  }

  if (renderingMode & OBJMESHRENDER_TEXTURE)
  {
    cgGLDisableProfile(FragmentProfile);
  }

  cgGLDisableProfile(VertexProfile);
}
*/

