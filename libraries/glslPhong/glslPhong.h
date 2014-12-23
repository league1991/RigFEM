/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "glslPhong" library , Copyright (C) 2014 USC                          *
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

/*
   Phong shading ("per-pixel lighting") in GLSL.
   Requires OpenGL 2.0 or higher
*/


#ifndef _GLSLPHONG_H_
#define _GLSLPHONG_H_

#include "openGL-headers.h"
#ifdef WIN32
#include "windows.h"
#endif

class GLSLPhong
{
public:
  GLSLPhong();
  ~GLSLPhong();

  void Enable(); // must be called each time one of the OpenGL lights is turned on/off
  void Disable();

protected:
  //GLint vertexShader;
  //GLint fragmentShader;
  //GLint program;
  GLint vertexShaders[256];
  GLint fragmentShaders[256];
  GLint programs[256];

  void checkError(GLint status, const char *msg);
  void checkShaderCompilation(GLint shaderID);

  static char vertexShaderStringAll [];
  static char vertexShaderStringPrologue [];
  static char vertexShaderStringCore [];
  static char vertexShaderStringEpilogue [];

  static char fragmentShaderStringAll [];
  static char fragmentShaderStringPrologue [];
  static char fragmentShaderStringCore [];
  static char fragmentShaderStringEpilogue [];
};

#endif

