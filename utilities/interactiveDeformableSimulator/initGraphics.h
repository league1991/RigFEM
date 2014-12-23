/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Interactive deformable object simulator" driver application,         *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin, Daniel Schroeder          *
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

#ifndef _INITGRAPHICS_H_
#define _INITGRAPHICS_H_

#ifdef WIN32
  #include <windows.h>
#endif

#include "openGL-headers.h"
#include "GL/glui.h"
#include "camera.h"

void initGLUT(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID);

void initCamera(double cameraRadius, 
				double cameraLongitude, double cameraLattitude,
				double focusPosX, double focusPosY, double focusPosZ,
				double camera2WorldScalingFactor,
				double * zNear, double * zFar, 
				SphericalCamera ** camera);

void initGraphics(int windowWidth, int windowHeight);

void setupLights();

// draw coordinate axes
void drawAxes(double axisLength);

void buildSphereDisplayList(GLuint * solidSphereList, GLuint * wireSphereList);

#endif

