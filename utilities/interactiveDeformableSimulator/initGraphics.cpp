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

#ifndef M_PI
  #define M_PI 3.141592653589793238462643
#endif

#include <math.h>
#include "initGraphics.h"

extern void displayFunction(void);
extern void idleFunction(void);
extern void reshape(int,int);
extern void handleMenu(int);
extern void keyboardFunction(unsigned char key, int x, int y);
extern void specialFunction(int key, int x, int y);
extern void mouseButtonActivityFunction(int button, int state, int x, int y);
extern void mouseMotionFunction(int x, int y);

// initialize GLUT
void initGLUT(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID)
{
  // Initialize GLUT.
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  //glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(windowWidth, windowHeight);
  *windowID = glutCreateWindow(windowTitle);

  // Setup GLUT callbacks.
  glutDisplayFunc(displayFunction); 

  GLUI_Master.set_glutIdleFunc(idleFunction);
  GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
  GLUI_Master.set_glutSpecialFunc(specialFunction);
  GLUI_Master.set_glutReshapeFunc(reshape);
  GLUI_Master.set_glutMouseFunc(mouseButtonActivityFunction);
    
  glutMotionFunc(mouseMotionFunction);

  // callback for mouse movement without any buttons pressed
  //glutPassiveMotionFunc(MouseNoDrag);
}

void initCamera(double cameraRadius, double cameraLongitude, double cameraLattitude, double focusPosX, double focusPosY, double focusPosZ, double camera2WorldScalingFactor, double * zNear, double * zFar, SphericalCamera ** camera)
{
  double focusPos[3] = {focusPosX, focusPosY, focusPosZ};

  *zNear = cameraRadius * 0.01;
  *zFar = cameraRadius * 100;

  double upPos[3] = {0,1,0};
  *camera = new SphericalCamera(cameraRadius, 1.0 * cameraLongitude / 360 * (2*M_PI), 1.0 * cameraLattitude / 360 * (2*M_PI), focusPos,  upPos, 0.05, camera2WorldScalingFactor); 

  (*camera)->SetOrigin(focusPos); 
}

void initGraphics(int windowWidth, int windowHeight)
{
   // clear to white
  glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to light blue
  //glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to gray
  //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);

  // clear to brown
  //glClearColor(255.0 / 256, 156.0 / 256, 17.0 / 256, 0.0);

  // clear to medical blue
  //glClearColor(148.0 / 256, 199.0 / 256, 211.0 / 256, 0.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_STENCIL_TEST);
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  reshape(windowWidth,windowHeight); 

  printf ("Graphics initialization complete.\n"); 
}

// draw the cartesian axes.
void drawAxes(double axisLength)
{
  glDisable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
    
  glBegin(GL_LINES);
  for (int i = 0; i<3; i++) 
  {
    float color[3] = { 0, 0, 0 };
    color[i] = 1.0;
    glColor3fv(color);
        
    float vertex[3] = {0, 0, 0};
    vertex[i] = axisLength;
    glVertex3fv(vertex);
    glVertex3f(0, 0, 0);
  } 
  glEnd();
}

