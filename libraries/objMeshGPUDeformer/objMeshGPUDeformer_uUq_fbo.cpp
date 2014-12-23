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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "glh_extensions.h"
#include "objMeshGPUDeformer_uUq_fbo.h"
//#include "GL/glext.h"

#if defined(linux)
  GLAPI void APIENTRY glGenFramebuffersEXT (GLsizei, GLuint *);
  GLAPI void APIENTRY glBindFramebufferEXT (GLenum, GLuint);
  GLAPI void APIENTRY glFramebufferTexture2DEXT (GLenum, GLenum, GLenum, GLuint, GLint);
  GLAPI GLenum APIENTRY glCheckFramebufferStatusEXT (GLenum);
#endif

#include "glh_extensions.h"

ObjMeshGPUDeformer_uUq_fbo::~ObjMeshGPUDeformer_uUq_fbo() {}

int ObjMeshGPUDeformer_uUq_fbo::InitRTT()
{
  if (InitExtensions() != 0)
  {
    printf("Init extensions failed.\n");
    return 1;
  }

  glGenFramebuffersEXT(1, &fbo);
  PrintGLerror("before binding texture to fbo");
  BindDynamicTextureToFBO();
  PrintGLerror("after binding texture to fbo");

  return 0;
}

void ObjMeshGPUDeformer_uUq_fbo::BindDynamicTextureToFBO()
{
  //printf("****Binding texture %d to fbo %d\n", dynamicTextureID, fbo);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
                            GL_TEXTURE_2D, vertexDeformationTextureID, 0);
  CheckFramebufferStatus();
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void ObjMeshGPUDeformer_uUq_fbo::EnableRTT()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, 
                            GL_TEXTURE_2D, vertexDeformationTextureID, 0);
}

void ObjMeshGPUDeformer_uUq_fbo::DisableRTT()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void ObjMeshGPUDeformer_uUq_fbo::BindRT()
{
}

void ObjMeshGPUDeformer_uUq_fbo::UnbindRT()
{
}

void ObjMeshGPUDeformer_uUq_fbo::CheckFramebufferStatus()
{
    GLenum status;
    status = (GLenum) glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch(status) {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
            printf("Framebuffer is complete.\n");
            break;
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
            printf("Unsupported framebuffer format\n");
            throw 51;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
            printf("Framebuffer incomplete, missing attachment\n");
            throw 51;
            break;
/*
        case GL_FRAMEBUFFER_INCOMPLETE_DUPLICATE_ATTACHMENT_EXT:
            printf("Framebuffer incomplete, duplicate attachment\n");
            throw 51;
            break;
*/
        case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
            printf("Framebuffer incomplete, attached images must have same dimensions\n");
            throw 51;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
            printf("Framebuffer incomplete, attached images must have same format\n");
            throw 51;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
            printf("Framebuffer incomplete, missing draw buffer\n");
            throw 51;
            break;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
            printf("Framebuffer incomplete, missing read buffer\n");
            throw 51;
            break;
        default:
            printf("Unknown framebuffer error.\n");
            throw 51;
    }
}

int ObjMeshGPUDeformer_uUq_fbo::InitExtensions()
{
  #ifdef WIN32
    if (!glh_init_extensions("GL_EXT_framebuffer_object "
                           "GL_ARB_multitexture "))
    {
      printf("Unable to load the following extension(s): %s\n", 
             glh_get_unsupported_extensions());
      printf("Try updating your graphics card driver.\n");
      return 1;
    }
    printf("Detected extensions: GL_EXT_framebuffer_object GL_ARB_multitexture\n");
  #endif

  return 0;
}

void ObjMeshGPUDeformer_uUq_fbo::SetDerivedData(void * data) 
{ 
  fbo = *(int*)data; 
  BindDynamicTextureToFBO();
}

