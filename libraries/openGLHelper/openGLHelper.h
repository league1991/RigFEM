/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "openGLHelper" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC   *
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
  Assorted OpenGL "helper" routines.
*/

#ifndef _OPENGLHELPER_H_
#define _OPENGLHELPER_H_

#ifdef WIN32
  #include <windows.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "macros.h"
#include "matrixIO.h"

#include "openGL-headers.h"

void OutputText(int x, int y, const char *string);

void UnitCube();
void UnitCubeWireframe();
void RenderWireframeBox(double bmin[3], double bmax[3]);

void DrawArrow( float px, float py, float pz,
    float nx, float ny, float nz,
    double arrowEndWidth, double arrowEndLength );

void DetermineCameraParameters(double centerX, double centerY, double centerZ, double modelRadius, double * focusX, double * focusY, double * focusZ, double * cameraRadius, double * zNear, double * zFar);

void JetColorMap(double x, double color[3]);

void RenderSphere(float x, float y, float z);
void BuildSphereDisplayList(GLuint * solidSphereList, GLuint * wireSphereList);
void TransparentSphere(GLuint solidSphereList, GLuint wireSphereList, double x, double y, double z, double radius);

char * DuplicateString(const char * s);

void PrintGLerror( const char *msg );

#endif

