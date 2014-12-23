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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ObjMeshGPUDeformer_uUq_pbuffer.h"

#ifdef WIN32

ObjMeshGPUDeformer_uUq_pbuffer::~ObjMeshGPUDeformer_uUq_pbuffer() {}

void ObjMeshGPUDeformer_uUq_pbuffer::EnableRTT()
{
  if( !wglMakeCurrent( pbuffer.hDC, pbuffer.hRC ) )
  {
    printf("Could not make the p-buffer's context current.");
    throw 1;
  }
}

void ObjMeshGPUDeformer_uUq_pbuffer::DisableRTT()
{
  // turn off pbuffer
  if( !wglMakeCurrent( hDC, hRC ) )
  {
    printf("Could not make the window's context current.");
    throw 2;
  }
}

void ObjMeshGPUDeformer_uUq_pbuffer::BindRT()
{
  if( !wglBindTexImageARB( pbuffer.hPBuffer, WGL_FRONT_LEFT_ARB ) )
  {
    printf("Could not bind p-buffer to render texture.");
    throw 3;
  }
}

void ObjMeshGPUDeformer_uUq_pbuffer::UnbindRT()
{
  // release the p-buffer from the dynamic "render-to" texture.
  if( !wglReleaseTexImageARB( pbuffer.hPBuffer, WGL_FRONT_LEFT_ARB ) )
  {
    printf("Could not release p-buffer from render texture.");
    throw 4;
  }
}

void ObjMeshGPUDeformer_uUq_pbuffer::DeleteRTT()
{
  wglDeleteContext(pbuffer.hRC);
  wglReleasePbufferDCARB( pbuffer.hPBuffer, pbuffer.hDC );
  wglDestroyPbufferARB( pbuffer.hPBuffer );
}

int ObjMeshGPUDeformer_uUq_pbuffer::InitRTT()
{
  // init WGL extensions

  wglGetExtensionsStringARB = (PFNWGLGETEXTENSIONSSTRINGARBPROC)wglGetProcAddress("wglGetExtensionsStringARB");
  char *ext = NULL;

  if( wglGetExtensionsStringARB )
    ext = (char*)wglGetExtensionsStringARB( wglGetCurrentDC() );
  else
  {
    printf("Unable to get address for wglGetExtensionsStringARB!");
    return 1;
  }

  //
  // WGL_ARB_pbuffer
  //

  if( strstr( ext, "WGL_ARB_pbuffer" ) == NULL )
  {
    printf("WGL_ARB_pbuffer extension was not found");
    return 1;
  }
  else
  {
    wglCreatePbufferARB    = (PFNWGLCREATEPBUFFERARBPROC)wglGetProcAddress("wglCreatePbufferARB");
    wglGetPbufferDCARB     = (PFNWGLGETPBUFFERDCARBPROC)wglGetProcAddress("wglGetPbufferDCARB");
    wglReleasePbufferDCARB = (PFNWGLRELEASEPBUFFERDCARBPROC)wglGetProcAddress("wglReleasePbufferDCARB");
    wglDestroyPbufferARB   = (PFNWGLDESTROYPBUFFERARBPROC)wglGetProcAddress("wglDestroyPbufferARB");
    wglQueryPbufferARB     = (PFNWGLQUERYPBUFFERARBPROC)wglGetProcAddress("wglQueryPbufferARB");

    if( !wglCreatePbufferARB || !wglGetPbufferDCARB || !wglReleasePbufferDCARB ||
      !wglDestroyPbufferARB || !wglQueryPbufferARB )
    {
      printf("One or more WGL_ARB_pbuffer functions were not found");
      return 1;
    }
  }

  //
  // WGL_ARB_pixel_format
  //

  if( strstr( ext, "WGL_ARB_pixel_format" ) == NULL )
  {
    printf("Error: WGL_ARB_pixel_format extension was not found");
    return 1;
  }
  else
  {
    wglChoosePixelFormatARB = (PFNWGLCHOOSEPIXELFORMATARBPROC)wglGetProcAddress("wglChoosePixelFormatARB");

    if( !wglChoosePixelFormatARB )
    {
      printf("Error: One or more WGL_ARB_pixel_format functions were not found");
      return 1;
    }
  }

  //
  // WGL_ARB_render_texture
  //

  if( strstr( ext, "WGL_ARB_render_texture" ) == NULL )
  {
    printf("WGL_ARB_render_texture extension was not found");
    return 1;
  }
  else
  {
    wglBindTexImageARB     = (PFNWGLBINDTEXIMAGEARBPROC)wglGetProcAddress("wglBindTexImageARB");
    wglReleaseTexImageARB  = (PFNWGLRELEASETEXIMAGEARBPROC)wglGetProcAddress("wglReleaseTexImageARB");
    wglSetPbufferAttribARB = (PFNWGLSETPBUFFERATTRIBARBPROC)wglGetProcAddress("wglSetPbufferAttribARB");

    if( !wglBindTexImageARB || !wglReleaseTexImageARB || !wglSetPbufferAttribARB )
    {
      printf("One or more WGL_ARB_render_texture functions were not found");
      return 1;
    }
  }

  hDC = wglGetCurrentDC();
  hRC = wglGetCurrentContext();

  //-------------------------------------------------------------------------
  // Create a p-buffer for off-screen rendering.
  //-------------------------------------------------------------------------
  int width = vertexDeformationTextureSize;
  int height = vertexDeformationTextureSize;

  pbuffer.hPBuffer = NULL;
  pbuffer.nWidth   = width;
  pbuffer.nHeight  = height;

  //
  // Define the minimum pixel format requirements we will need for our 
  // p-buffer. A p-buffer is just like a frame buffer, it can have a depth 
  // buffer associated with it and it can be double buffered.
  //

  int pf_attr[] =
  {
    WGL_SUPPORT_OPENGL_ARB, TRUE,       // P-buffer will be used with OpenGL
      WGL_DRAW_TO_PBUFFER_ARB, TRUE,      // Enable render to p-buffer
      WGL_BIND_TO_TEXTURE_RGBA_ARB, TRUE, // P-buffer will be used as a texture
      WGL_RED_BITS_ARB, 32,                // At least 8 bits for RED channel
      WGL_GREEN_BITS_ARB, 32,              // At least 8 bits for GREEN channel
      WGL_BLUE_BITS_ARB, 32,               // At least 8 bits for BLUE channel
      WGL_ALPHA_BITS_ARB, 32,              // At least 8 bits for ALPHA channel
      WGL_DEPTH_BITS_ARB, 0,             // At least 0 bits for depth buffer
      WGL_PIXEL_TYPE_ARB, WGL_TYPE_RGBA_FLOAT_ATI,
      WGL_DOUBLE_BUFFER_ARB, FALSE,       // We don't require double buffering
      0                                   // Zero terminates the list
  };

  unsigned int count = 0;
  int pixelFormat;
  wglChoosePixelFormatARB( hDC, (const int*)pf_attr, NULL, 1, &pixelFormat, &count);

  if( count == 0 )
  {
    printf("Could not find an acceptable pixel format!");
    return 1;
  }

  //
  // Set some p-buffer attributes so that we can use this p-buffer as a
  // 2D RGBA texture target.
  //

  int pb_attr[] =
  {
    WGL_TEXTURE_FORMAT_ARB, WGL_TEXTURE_RGBA_ARB, // Our p-buffer will have a texture format of RGBA
      WGL_TEXTURE_TARGET_ARB, WGL_TEXTURE_2D_ARB,   // Of texture target will be GL_TEXTURE_2D
      0                                             // Zero terminates the list
  };

  //
  // Create the p-buffer...
  //

  pbuffer.hPBuffer = wglCreatePbufferARB( hDC, pixelFormat, 
    pbuffer.nWidth, pbuffer.nHeight, pb_attr );
  pbuffer.hDC      = wglGetPbufferDCARB( pbuffer.hPBuffer );
  pbuffer.hRC      = wglCreateContext( pbuffer.hDC );

  if( !pbuffer.hPBuffer )
  {
    printf("Error: could not create the p-buffer");
    return 1;
  }

  int queryHeight;
  int queryWidth;
  wglQueryPbufferARB( pbuffer.hPBuffer, WGL_PBUFFER_WIDTH_ARB, &queryHeight );
  wglQueryPbufferARB( pbuffer.hPBuffer, WGL_PBUFFER_WIDTH_ARB, &queryWidth );

  if( queryHeight != pbuffer.nHeight || queryWidth != pbuffer.nWidth )
  {
    printf("The width and height of the created p-buffer don't match the requirements!");
    return 1;
  }

  return 0;
}

// initialize static members
HDC ObjMeshGPUDeformer_uUq_pbuffer::hDC = NULL;
HGLRC ObjMeshGPUDeformer_uUq_pbuffer::hRC = NULL;
PFNWGLGETEXTENSIONSSTRINGARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglGetExtensionsStringARB = NULL;
PFNWGLCREATEPBUFFERARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglCreatePbufferARB = NULL;
PFNWGLGETPBUFFERDCARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglGetPbufferDCARB = NULL;
PFNWGLRELEASEPBUFFERDCARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglReleasePbufferDCARB = NULL;
PFNWGLDESTROYPBUFFERARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglDestroyPbufferARB = NULL;
PFNWGLQUERYPBUFFERARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglQueryPbufferARB = NULL;
PFNWGLCHOOSEPIXELFORMATARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglChoosePixelFormatARB = NULL;
PFNWGLBINDTEXIMAGEARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglBindTexImageARB = NULL;
PFNWGLRELEASETEXIMAGEARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglReleaseTexImageARB = NULL;
PFNWGLSETPBUFFERATTRIBARBPROC ObjMeshGPUDeformer_uUq_pbuffer::wglSetPbufferAttribARB = NULL;

#endif

