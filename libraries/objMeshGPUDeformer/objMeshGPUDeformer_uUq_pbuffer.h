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

// this class works under Windows only

#ifndef _OBJMESHGPUDEFORMER_UUQ_PBUFFER_H_
#define _OBJMESHGPUDEFORMER_UUQ_PBUFFER_H_

#ifdef WIN32

#ifdef _MSC_VER
  #pragma comment( lib, "cg.lib" )
  #pragma comment( lib, "cgGL.lib" )
#endif

//#include <Windows.h>
#include <gl/glew.h>
#include <gl/GLEXT.H>
#include <gl/wglew.h>
#include "objMeshGPUDeformer_uUq.h"
//#include <GL/wglext.h>

class ObjMeshGPUDeformer_uUq_pbuffer : public ObjMeshGPUDeformer_uUq
{
public:
  virtual ~ObjMeshGPUDeformer_uUq_pbuffer();

protected:

  virtual void EnableRTT();
  virtual void DisableRTT();
  virtual void BindRT();
  virtual void UnbindRT();
  virtual int InitRTT();
  virtual void DeleteRTT();

  virtual void * GetDerivedData() { return &pbuffer; }
  virtual void SetDerivedData(void * data) { pbuffer = *(PBUFFER*)data; }

  // pbuffer and wgl data
  typedef struct   
  {
    HPBUFFERARB hPBuffer;
    HDC         hDC;
    HGLRC       hRC;
    int         nWidth;
    int         nHeight;
  } PBUFFER;

  PBUFFER pbuffer;

  static HDC	hDC;
  static HGLRC  hRC;

  // WGL_ARB_extensions_string
  static PFNWGLGETEXTENSIONSSTRINGARBPROC wglGetExtensionsStringARB;

  // WGL_ARB_pbuffer
  static PFNWGLCREATEPBUFFERARBPROC    wglCreatePbufferARB;
  static PFNWGLGETPBUFFERDCARBPROC     wglGetPbufferDCARB;
  static PFNWGLRELEASEPBUFFERDCARBPROC wglReleasePbufferDCARB;
  static PFNWGLDESTROYPBUFFERARBPROC   wglDestroyPbufferARB;
  static PFNWGLQUERYPBUFFERARBPROC     wglQueryPbufferARB;

  // WGL_ARB_pixel_format
  static PFNWGLCHOOSEPIXELFORMATARBPROC wglChoosePixelFormatARB;

  // WGL_ARB_render_texture
  static PFNWGLBINDTEXIMAGEARBPROC     wglBindTexImageARB;
  static PFNWGLRELEASETEXIMAGEARBPROC  wglReleaseTexImageARB;
  static PFNWGLSETPBUFFERATTRIBARBPROC wglSetPbufferAttribARB;
};

#endif

#endif
