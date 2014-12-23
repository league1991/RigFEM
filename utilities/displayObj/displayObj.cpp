/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * OBJ mesh visualization utility                                        *
 * Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                            *
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
  OBJ visualization utility.
  Renders an obj file (or multiple objs) to the screen using OpenGL.
  A group can be selected/highlighted (ctrl + LMB);
*/

#ifdef WIN32
  #pragma warning(disable : 4996)
  #pragma warning(disable : 4267)
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <set>

#ifdef WIN32
#include "GL\glew.h"
#include <windows.h>
#ifndef M_PI
#define M_PI 3.141592654 
#endif
#endif

#include "openGL-headers.h"
#include "objMeshRender.h"
#include "getopts.h"
#include "matrixIO.h"
#include "GL/glui.h"
#include "performanceCounter.h"
#include "openGLHelper.h"
#include "lighting.h"
#include "glslPhong.h"
#include "camera.h"
#include "loadList.h"

std::set<int> * selectedVerticesSet;  // each object will have a set storing selected vertices
std::set<int> * selectedGroupsSet; // each object will have a set storing selected groups

int currentSelectedObject = -1;
int currentSelectedGroup = -1;  // on current selected object
int currentSelectedVertex = -1; // on current selected object

int sprite = 0;
int windowWidth, windowHeight;
SphericalCamera * camera;

std::map<unsigned int,ObjMesh *> objMeshes;
std::map<unsigned int,ObjMeshRender *> objMeshRenders;
std::map<unsigned int,GLuint> displayLists;
std::map<unsigned int,GLuint> displayListsWithEdges;

char * objMeshname;
char fpsString[32] = "--";
double radius;
Vec3d lightBoxMin, lightBoxMax;
int numLights =2;
int showSpecularLight = 0;
int displayRuntimeInfo = 0;

char screenShotFilename[FILENAME_MAX] = "__default";
#define SCREEN_SHOT_DEFAULT_WIDTH 800
#define SCREEN_SHOT_DEFAULT_HEIGHT 800
bool windowResizedFlag = false;

bool autoRotateModelFlag = false; // rotate model so the longest edge of the bounding box aligns the Y axis
bool rotateZUpToYUpFlag = false;

int renderObject = 1;
int renderLights = 1;
int enableGLSLPhong = 1;
int showAxes = 1;
int lightOn[8] = {1, 0, 0, 0, 0, 1, 0, 1 };
int showNormals = 0;
float lightIntensity = 0.5;
double lastPrintOuttime = 0.0;
char lightingConfigFilename[4096] = "__none";
Lighting * lighting = NULL;

PerformanceCounter performanceCounter;
double lastTime = 0;

bool useFieldColors;
bool renderMaterials;
std::vector<Vec3d> fieldColors;
int numTextures;
int renderTextures = 1;

int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;
int g_iAlt = 0;
int g_iShift = 0;
int g_iCtrl = 0;
int mouseUpdates=0;
float renderNormalLength = 0.05;

bool renderExternalVertices;
int numExternalVertices;
double * externalVertices;

GLUI * glui;

GLSLPhong * glslPhong = NULL;

int enableVertexSelection = 1;

int winID;

bool displayEdges = false;

void PrintSelectedVertices()
{
  for(unsigned int objectIndex=0; objectIndex<objMeshes.size(); objectIndex++)
    if (selectedVerticesSet[objectIndex].size() != 0)
    {
      printf("Number of selected vertices on object %d: %d | 0-indexed:\n", (int)objectIndex, (int)selectedVerticesSet[objectIndex].size());    
      for(std::set<int>::iterator iter = selectedVerticesSet[objectIndex].begin(); iter != selectedVerticesSet[objectIndex].end(); iter++)
        printf("%d ", *iter);
      printf("\n");
      printf("Number of selected vertices on object %d: %d | 1-indexed:\n", (int)objectIndex, (int)selectedVerticesSet[objectIndex].size());    
      for(std::set<int>::iterator iter = selectedVerticesSet[objectIndex].begin(); iter != selectedVerticesSet[objectIndex].end(); iter++)
        printf("%d ", *iter + 1);    
      printf("\n");
    }
}

void MakeDisplayList(unsigned int fileIndex)
{
  int renderMode;
  if (renderMaterials)
  {
    renderMode = OBJMESHRENDER_SMOOTH | OBJMESHRENDER_MATERIAL;
    if ((numTextures > 0) && (renderTextures))
    {
      renderMode |= OBJMESHRENDER_TEXTURE;
      if (objMeshRenders[fileIndex]->maxBytesPerPixelInTextures() == 4)
        renderMode |= OBJMESHRENDER_TRANSPARENCY;
    }
  }
  else
    renderMode = OBJMESHRENDER_CUSTOMCOLOR;

  int geometryMode = OBJMESHRENDER_TRIANGLES;
  GLuint objMeshList = objMeshRenders[fileIndex]->createDisplayList(geometryMode, renderMode);
  displayLists.insert(std::make_pair(fileIndex, objMeshList));

  geometryMode = OBJMESHRENDER_EDGES;
  renderMode = OBJMESHRENDER_NONE;
  GLuint objMeshListWithEdges = objMeshRenders[fileIndex]->createDisplayList(geometryMode, renderMode);
  displayListsWithEdges.insert(std::make_pair(fileIndex, objMeshListWithEdges));
}

void MakeDisplayLists()
{
  for(unsigned int i=0; i < objMeshes.size(); i++)
    MakeDisplayList(i);
}

void EraseDisplayLists()
{
  std::map<unsigned int,GLuint> :: iterator iter;

  for(iter=displayLists.begin(); iter != displayLists.end(); iter++)
    glDeleteLists(displayLists[iter->second],1);

  for(iter=displayListsWithEdges.begin(); iter != displayListsWithEdges.end(); iter++)
    glDeleteLists(displayListsWithEdges[iter->second],1);
}

void SetupLights()
{
  if (lighting != NULL)
    lighting->LightScene();
  else
  {
    GLfloat aGa[] = { lightIntensity, lightIntensity, lightIntensity, 1.0 };
  
    GLfloat lKa[] = { lightIntensity, lightIntensity, lightIntensity, 1.0 };
    GLfloat lKd[] = { lightIntensity, lightIntensity, lightIntensity, 1.0 };
    GLfloat lKs[] = { lightIntensity, lightIntensity, lightIntensity, 1.0 };
  
    if (showSpecularLight == 0) 
    {
      lKs[0] = 0.0; lKs[1] = 0.0; lKs[2] = 0.0;
    }
  
    // light positions and directions 
    //GLfloat lP0[4] = { 0.2, 0.2, -0.2, 1.0 }; 
    GLfloat lP0[] = { lightBoxMin[0], lightBoxMin[1], lightBoxMin[2], 1.0 };
    GLfloat lP1[] = { lightBoxMax[0], lightBoxMin[1], lightBoxMin[2], 1.0 };
    GLfloat lP2[] = { lightBoxMax[0], lightBoxMax[1], lightBoxMin[2], 1.0 };
    GLfloat lP3[] = { lightBoxMin[0], lightBoxMax[1], lightBoxMin[2], 1.0 };
    GLfloat lP4[] = { lightBoxMin[0], lightBoxMin[1], lightBoxMax[2], 1.0 };
    GLfloat lP5[] = { lightBoxMax[0], lightBoxMin[1], lightBoxMax[2], 1.0 };
    GLfloat lP6[] = { lightBoxMax[0], lightBoxMax[1], lightBoxMax[2], 1.0 };
    GLfloat lP7[] = { lightBoxMin[0], lightBoxMax[1], lightBoxMax[2], 1.0 };
  
    // set up 
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
    glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    //glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  
    #define LIGHTSETUP(i)\
    glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
    glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa);\
    glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd);\
    glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs);\
  
    LIGHTSETUP(0);
    LIGHTSETUP(1);
    LIGHTSETUP(2);
    LIGHTSETUP(3);
    LIGHTSETUP(4);
    LIGHTSETUP(5);
    LIGHTSETUP(6);
    LIGHTSETUP(7);
  
    #define TURNONLIGHT(i)\
     if(lightOn[i])\
       glEnable(GL_LIGHT##i);\
     else\
       glDisable(GL_LIGHT##i);\
    
    TURNONLIGHT(0);
    TURNONLIGHT(1);
    TURNONLIGHT(2);
    TURNONLIGHT(3);
    TURNONLIGHT(4);
    TURNONLIGHT(5);
    TURNONLIGHT(6);
    TURNONLIGHT(7);
  
    glEnable(GL_LIGHTING);
  }
}
  
void DoIdle()
{
  glutSetWindow(winID);

  performanceCounter.StopCounter();
  double time = performanceCounter.GetElapsedTime();

  double fps = 1.0 / (time - lastTime);
  lastTime = time;

  if (time - lastPrintOuttime > 1.0)
  {
    lastPrintOuttime = time;
    sprintf(fpsString,"%.1f Hz",fps);
    char titleString[4096];
    sprintf(titleString, "%s | %s", objMeshname, fpsString);
    glutSetWindowTitle(titleString);
  }

/*
  // save screen shot
  if (strcmp(screenShotFilename, "__default") != 0)
  {
    if (sprite >= 3) // wait for three frames (safety issue)
    {
       char filenameExtension[FILENAME_MAX];
       strcpy(filenameExtension, &screenShotFilename[strlen(screenShotFilename) - 3]);
       for(int i=0; i<(int)strlen(filenameExtension); i++)
         filenameExtension[i] = tolower(filenameExtension[i]);

      if (strcmp(filenameExtension, "jpg") == 0)
        SaveScreenshot(screenShotFilename, SCREEN_SHOT_DEFAULT_WIDTH, SCREEN_SHOT_DEFAULT_HEIGHT);

      if (strcmp(filenameExtension, "tif") == 0)
        SaveScreenshotTiff(screenShotFilename, SCREEN_SHOT_DEFAULT_WIDTH, SCREEN_SHOT_DEFAULT_HEIGHT);    
 
      if (strcmp(filenameExtension, "ppm") == 0)
        SaveScreenshotPPM(screenShotFilename, SCREEN_SHOT_DEFAULT_WIDTH, SCREEN_SHOT_DEFAULT_HEIGHT);
      printf("-------------------------------------------------------\n");
    }
    if (sprite >= 5)  // screenshot will be saved three times at this point, however, if we don't do this, we sometimes get a totally black picture!! 
     exit(0);
  }
*/

  sprite++;
  glutPostRedisplay();
}

void Reshape(int x, int y)
{
  glViewport(0,0,x,y);

  glMatrixMode(GL_PROJECTION); // Select The Projection Matrix
  glLoadIdentity(); // Reset The Projection Matrix

  windowWidth = x;
  windowHeight = y;

  // Calculate The Aspect Ratio Of The Window
  gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, 0.002 * radius, 10 * radius);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void InitGLUI()
{
  glui = GLUI_Master.create_glui("Lights",0,650,0);

  if (lighting == NULL)
  {
    for(int i=0; i<8; i++)
    {
      char s[20];
      sprintf(s,"Light %d",i);
      glui->add_checkbox( s, &lightOn[i] );
    }

    GLUI_Spinner * lightIntensitySpinner = glui->add_spinner("Light intensity:", GLUI_SPINNER_FLOAT, &lightIntensity);
    lightIntensitySpinner->set_float_limits(0,1E6);
    lightIntensitySpinner->set_speed(0.2);
    glui->add_checkbox("Specular term", &showSpecularLight);
  }

  glui->add_checkbox("Render normals", &showNormals);
  GLUI_Spinner * normalLengthSpinner = glui->add_spinner("Normal length:", GLUI_SPINNER_FLOAT, &renderNormalLength);
  normalLengthSpinner->set_float_limits(0,1E6);
  normalLengthSpinner->set_speed(0.2);
}

void Initialize()
{
  sprite = 0;

  // rotate model if needed
  if (autoRotateModelFlag)
  {
    // compute the bounding box of the input obj mesh
    for(int objectIndex=0; objectIndex<(int)objMeshes.size(); objectIndex++)
    {
      Vec3d bbmin, bbmax;
      objMeshes[objectIndex]->getBoundingBox(1.0, &bbmin, &bbmax);
      Vec3d edge = bbmax - bbmin;

      if ((edge[0] > edge[1]) && (edge[0] > edge[2]))  // edge on x direction is the longest
      {
        Vec3d translation = Vec3d(0, 0, 0);
        Mat3d rotation(0, -1, 0, 1, 0, 0, 0, 0, 1);  // rotate so that, x_new = -y; y_new = x; z_new = z;
        objMeshes[objectIndex]->transformRigidly(translation, rotation);
      }

      if ((edge[2] > edge[0]) && (edge[2] > edge[1]))  // edge on z direction is the longest
      {
        Vec3d translation = Vec3d(0, 0, 0);
        Mat3d rotation(1, 0, 0, 0, 0, 1, 0, -1, 0); // rotate so that, x_new = x; y_new = z; z_new = -y;
        objMeshes[objectIndex]->transformRigidly(translation, rotation);
      }
    }  // for
  }  // if 

  if (rotateZUpToYUpFlag)
  {
    for(int objectIndex=0; objectIndex<(int)objMeshes.size(); objectIndex++)
    {
      Vec3d translation = Vec3d(0, 0, 0);
      Mat3d rotation(1, 0, 0, 0, 0, 1, 0, -1, 0); // rotate so that, x_new = x; y_new = z; z_new = -y;
      objMeshes[objectIndex]->transformRigidly(translation, rotation);
    }  // for
  }

  // compute bounding box of the scene
  Vec3d sceneBoxMin;
  Vec3d sceneBoxMax;
  objMeshes[0]->getBoundingBox(1.0, &sceneBoxMin, &sceneBoxMax);
  for(int objectIndex=1; objectIndex<(int)objMeshes.size(); objectIndex++)
  {
    Vec3d bmin, bmax;
    (objMeshes[objectIndex])->getCubicBoundingBox(1.25, &bmin, &bmax); 
    for(int dim=0; dim<3; dim++)
    {
      if (sceneBoxMin[dim] > bmin[dim]) sceneBoxMin[dim] = bmin[dim];
      if (sceneBoxMax[dim] < bmax[dim]) sceneBoxMax[dim] = bmax[dim];
    }
  }

  Vec3d center;
  center = 0.5 * (sceneBoxMin + sceneBoxMax);

  radius = 2.5 * len(sceneBoxMax-center);

  double centerV[3] = {center[0], center[1], center[2]};
  double upVector[3] = {0,1,0};

  if (strcmp(lightingConfigFilename, "__none") != 0)
    lighting = new Lighting(lightingConfigFilename);
  else
    lighting = NULL;

  // init camera for 3d viewing
  // bogusFlag, r, Theta, Phi, focusX, focusY, focusZ, upX, upY, upZ, sensitivity, camera2worldScalingFactor
  //camera = new SphericalCamera(radius, 1.0 * 270 / 360 * (2*PI), 1.0 * 30 / 360 * (2*PI), centerV, upVector,  0.05, 1);
  camera = new SphericalCamera(radius, 1.0 * 270 / 360 * (2*PI), 1.0 * 10 / 360 * (2*PI), centerV, upVector,  0.05, 1);

  // clear to white
  glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to light blue
  //glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to gray
  //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_STENCIL_TEST);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  //glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);

  //glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_NORMAL_ARRAY);

  //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Really Nice Perspective Calculations

  GLint AAbufs, AAsamples;
  glGetIntegerv(GL_SAMPLE_BUFFERS, &AAbufs);
  glGetIntegerv(GL_SAMPLES, &AAsamples);
  printf("AA: bufs=%d samples=%d\n", AAbufs, AAsamples);
  glEnable(GL_MULTISAMPLE);

  MakeDisplayLists(); // cache the geometry

  Reshape(windowWidth,windowHeight); // calls gluPerspective, among other things

  glslPhong = new GLSLPhong();

  if (strcmp(screenShotFilename, "__default") != 0)
  {
    showAxes = 0;
    showNormals = 0;
    renderLights = 0;
    displayRuntimeInfo = 0;
  }  

  return;
}

void MouseNoDrag(int x, int y)
{
  mouseUpdates++;

/*
  int mouseDeltaX = x-g_vMousePos[0];
  int mouseDeltaY = y-g_vMousePos[1];

  if (mouseUpdates == 1) 
  {
    mouseDeltaX = 0;
    mouseDeltaY = 0;
  }
*/

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void MouseDrag(int x, int y)
{ 
  int mouseDeltaX = x-g_vMousePos[0];
  int mouseDeltaY = y-g_vMousePos[1];

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
    
  // we moved the camera...
  if (g_iRightMouseButton) // right mouse button handle camera rotations
  {
    camera->MoveRight(mouseDeltaX);
    camera->MoveUp(mouseDeltaY);
  }

  if ((g_iMiddleMouseButton) || (g_iAlt && g_iLeftMouseButton)) // handle camera zoom
    camera->ZoomIn(mouseDeltaY*radius*0.05);

/*
  if (g_iMiddleMouseButton)
  {
    camera->MoveIn(mouseDeltaY*radius*0.025);
    camera->ZoomIn(mouseDeltaY*radius*0.025);
  }
*/
}

void MouseButtonActivity(int button, int state, int x, int y)
{
  g_iShift = glutGetModifiers() & GLUT_ACTIVE_SHIFT;
  g_iAlt = glutGetModifiers() & GLUT_ACTIVE_ALT;
  g_iCtrl = glutGetModifiers() & GLUT_ACTIVE_CTRL;

  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      if (enableVertexSelection && (g_iLeftMouseButton))
      {
        GLdouble model[16];
        glGetDoublev (GL_MODELVIEW_MATRIX, model);

        GLdouble proj[16];
        glGetDoublev (GL_PROJECTION_MATRIX, proj);

        GLint view[4];
        glGetIntegerv (GL_VIEWPORT, view);

        int winX = x;
        int winY = view[3]-1-y;

        float zValue;
        glReadPixels(winX,winY,1,1, GL_DEPTH_COMPONENT, GL_FLOAT, &zValue);

        GLubyte stencilValue;
        glReadPixels(winX, winY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, &stencilValue);

        GLdouble worldX, worldY, worldZ;
        gluUnProject (winX, winY, zValue, model, proj, view, &worldX, &worldY, &worldZ);

        currentSelectedObject = -1;
        currentSelectedGroup = -1;
        currentSelectedVertex = -1;
        if(g_iShift)  // vertex selection mode
        {
          currentSelectedObject = stencilValue - 1;
          if (stencilValue != 0)
          {
            double distance;
            int clickedVertex = objMeshes[currentSelectedObject]->getClosestVertex(Vec3d(worldX, worldY, worldZ), &distance);
            printf("Clicked on vertex %d of object %d\n", clickedVertex, currentSelectedObject);
            if (selectedVerticesSet[currentSelectedObject].find(clickedVertex) != selectedVerticesSet[currentSelectedObject].end())
              selectedVerticesSet[currentSelectedObject].erase(clickedVertex);
            else
              selectedVerticesSet[currentSelectedObject].insert(clickedVertex);
            currentSelectedVertex = clickedVertex;
          }
          else
          {
            printf("Clicked on empty stencil: %d.\n", stencilValue);
          }
          PrintSelectedVertices();
        }

        if (g_iCtrl)  // group selection mode
        {
          currentSelectedObject = stencilValue - 1;
          if (stencilValue != 0)
          {
            int clickedVertex = objMeshes[currentSelectedObject]->getClosestVertex(Vec3d(worldX, worldY, worldZ));
            // find which group this vertex is belongs to:
            bool vertexFoundFlag = false;
            for(int groupIndex=0; !vertexFoundFlag && groupIndex<(int)objMeshes[currentSelectedObject]->getNumGroups(); groupIndex++)
            {
              const ObjMesh::Group * groupHandle = objMeshes[currentSelectedObject]->getGroupHandle(groupIndex);
              for(int faceIndex=0; !vertexFoundFlag && faceIndex<(int)groupHandle->getNumFaces(); faceIndex++)
              {
                const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(faceIndex);
                for(int vertexIndex=0; !vertexFoundFlag && vertexIndex<(int)faceHandle->getNumVertices(); vertexIndex++)
                  if ((int)(faceHandle->getVertexHandle(vertexIndex)->getPositionIndex()) == clickedVertex)
                  {
                    vertexFoundFlag = true;
                    currentSelectedGroup = groupIndex;
                    printf("Clicked on group %d of object %d\n", currentSelectedGroup, currentSelectedObject);
                    //if (selectedGroupsSet[currentSelectedObject].find(currentSelectedGroup) != selectedGroupsSet[currentSelectedObject].end())
                    //  selectedGroupsSet[currentSelectedObject].erase(currentSelectedGroup);
                    //else
                    //  selectedGroupsSet[currentSelectedObject].insert(currentSelectedGroup);
                  }
              }  // for faceIndex
            }  // for groupIndex
          }  // if (stencilValue != 0)    
          else
          {
            currentSelectedGroup = -1;
            printf("Clicked on empty stencil: %d.\n", stencilValue);
          }         

          if (currentSelectedGroup >= 0)
            displayRuntimeInfo = 1;

        }  // if (g_iCtrl)
      }

      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void Recycle(void)
{
  delete [] selectedVerticesSet;
  delete [] selectedGroupsSet;
}
void KeyboardFunc (unsigned char key, int x, int y)
{ 
  double cameraX,cameraY,cameraZ;
  const double earSeparation = 0.2;
  double leftEarX,leftEarY,leftEarZ, rightEarX,rightEarY,rightEarZ;
  GLfloat m[16];

  switch (key)
  {
    case 27:
      EraseDisplayLists();
      glutDestroyWindow(winID);
      Recycle();
      exit(0);
      break;

    case 'w':
      displayEdges = !displayEdges;
      break;

    case 'e':
      renderObject = 1 - renderObject;
      break;

    case 'n':
      showNormals = 1 - showNormals;
      if(renderMaterials)
        glui->sync_live();
      break;

    case 'f':
      displayRuntimeInfo = !displayRuntimeInfo;
    break;

    case '\\': 
      camera->Reset();
      break;

    case ']':
    {
      char pos[96] = "cameraPos.txt";
      camera->SavePosition(pos);
    }
    break;

    case '[':
    {
      char pos[96] = "cameraPos.txt";
      camera->LoadPosition(pos);
    }
    break;

    case 'a': 
      showAxes = 1 - showAxes;
      break;

    case 't': 
      renderTextures = 1 - renderTextures;
      break;

    case 'L': 
      renderLights = 1 - renderLights;
      break;

    case 'P': 
      enableGLSLPhong = 1 - enableGLSLPhong;
      printf("GLSL Phong enabled: %d\n", enableGLSLPhong);
      break;

    case 'i':
      camera->GetAbsWorldPosition(cameraX,cameraY,cameraZ);
      printf("Camera is positioned at: %G %G %G\n",cameraX,cameraY,cameraZ);
      printf("Camera radius is: %G \n",camera->GetRadius());
      printf("Camera Phi is: %G \n",180.0/M_PI*camera->GetPhi());
      printf("Camera Theta is: %G \n",180.0/M_PI*camera->GetTheta());
      camera->GetStereoPosition(earSeparation,&leftEarX,&leftEarY,&leftEarZ,
        &rightEarX,&rightEarY,&rightEarZ);
      printf("Stereo: left ear is at: %G %G %G\n",leftEarX,leftEarY,leftEarZ);
      printf("Stereo: right ear is at: %G %G %G\n",rightEarX,rightEarY,rightEarZ);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
        glLoadIdentity();
        camera->Look(); // calls gluLookAt
        glGetFloatv (GL_MODELVIEW_MATRIX, m);
      glPopMatrix();
      printf("Modelview matrix is:\n");
      printf("%f %f %f %f\n",m[0],m[4],m[8],m[12]);
      printf("%f %f %f %f\n",m[1],m[5],m[9],m[13]);
      printf("%f %f %f %f\n",m[2],m[6],m[10],m[14]);
      printf("%f %f %f %f\n",m[3],m[7],m[11],m[15]);
      break;


    #define LIGHTCASE(i)\
      case 48+i:\
      lightOn[i] = !lightOn[i];\
      break;

    LIGHTCASE(0);
    LIGHTCASE(1);
    LIGHTCASE(2);
    LIGHTCASE(3);
    LIGHTCASE(4);
    LIGHTCASE(5);
    LIGHTCASE(6);
    LIGHTCASE(7);
  }
}

void SpecialKeysFunc (int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      camera->MoveFocusRight(+0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_RIGHT:
      camera->MoveFocusRight(-0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_DOWN:
      camera->MoveFocusUp(+0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_UP:
      camera->MoveFocusUp(-0.1 * camera->GetRadius());
    break;

    case GLUT_KEY_PAGE_UP:
    break;

    case GLUT_KEY_PAGE_DOWN:
    break;

    case GLUT_KEY_HOME:
    break;

    case GLUT_KEY_END:
    break;

    case GLUT_KEY_INSERT:
    break;
/*
    case GLUT_KEY_DELETE:
    break;
*/

    default:
      break;
  }
}

//font is e.g. GLUT_BITMAP_9_BY_15
void print_bitmap_string(float x,float y, float z, void *font, char* s)
{
  glRasterPos3f(x,y,z);
  if (s && strlen(s)) 
  {
    while (*s) 
    {
      glutBitmapCharacter(font, *s);
      s++;
    }
  }
}

void DisplayRuntimeInformation(void)
{ 
  // 2D rasterized information
  // bitmap routines
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  if (displayRuntimeInfo)
  {
    glColor3f(0.0f, 0.0f, 0.0f);
    gluOrtho2D (0, windowWidth, 0, windowHeight);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    char s[4096];
    int lineCount = 1;
    int X1 = 10;
    int Y1;
    //int X1 = -1 + 2.0 * x1 / windowWidth;
    //int Y1 = -1 + 2.0 * y1 / windowHeight;
    glColor3f(0,0,0);
    s[0] = '\0';
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    if (currentSelectedObject >= 0)
      sprintf(s, "Selected Object: %d", currentSelectedObject);
    else
      sprintf(s, "Selected Object: unknown");
    Y1 = windowHeight - 20 * lineCount;
    print_bitmap_string(X1,Y1,-0.9,GLUT_BITMAP_9_BY_15,s);
    lineCount++;

    if (currentSelectedGroup >= 0 && currentSelectedObject >= 0)
      sprintf(s, "Selected Group: %d", currentSelectedGroup);
    else
      sprintf(s, "Selected Group: unknown");
    Y1 = windowHeight - 20 * lineCount;
    print_bitmap_string(X1,Y1,-0.9,GLUT_BITMAP_9_BY_15,s);
    lineCount++;

    if (currentSelectedVertex >= 0 && currentSelectedObject >= 0)
    {
      Vec3d pos = objMeshes[currentSelectedObject]->getPosition(currentSelectedVertex);
      sprintf(s, "Selected vertex pos: %.02f, %.02f, %.02f", pos[0], pos[1], pos[2]);
    }
    else
      sprintf(s, "Selected vertex pos: unknown");
    Y1 = windowHeight - 20 * lineCount;
    print_bitmap_string(X1,Y1,-0.9,GLUT_BITMAP_9_BY_15,s);
    lineCount++;

  }  // if (!displayHelperInfo)
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  // end of bitmap
}

void Display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  // Setup model transformations
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();

  camera->Look(); // calls gluLookAt

  glDisable(GL_LIGHTING);

  // turns off stencil modifications
  glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

  glLineWidth(2.0);
  if (showAxes)
  {
    glBegin(GL_LINES);
    glColor3f(1,0,0);
      glVertex3f(0,0,0);
      glVertex3f(1,0,0);
    glColor3f(0,1,0);
      glVertex3f(0,0,0);
      glVertex3f(0,1,0);
    glColor3f(0,0,1);
      glVertex3f(0,0,0);
      glVertex3f(0,0,1);
    glEnd();
  }

  if (showNormals)
  {
    glColor3f(1,0,0);
    for(unsigned int i=0; i<objMeshes.size(); i++)
      objMeshRenders[i]->renderNormals(renderNormalLength);
  }

  if (renderMaterials)
    SetupLights();

  // turn on stencil

//  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
//  glStencilFunc(GL_ALWAYS, 1, ~(0u));

  glStencilOp(GL_KEEP, GL_REPLACE, GL_REPLACE);
  glStencilFunc(GL_ALWAYS, 1, ~(0u));

  if (enableGLSLPhong)
    glslPhong->Enable();

  if (renderObject == 1)
  {
    for(unsigned int i=0; i<objMeshes.size(); i++)
    {
      glStencilFunc(GL_ALWAYS, i+1, ~(0u));
      unsigned int objMeshList = displayLists[i];
      glEnable(GL_LIGHTING);
      glCallList(objMeshList); 
    }
  }

  if (displayEdges)
  {
    for(unsigned int i=0; i<objMeshes.size(); i++)
    {
      glStencilFunc(GL_ALWAYS, i+1, ~(0u));
      unsigned int objMeshListWithEdges = displayListsWithEdges[i];
      glDisable(GL_LIGHTING);
      glColor3f(0,0,0);
      glslPhong->Disable();
      glCallList(objMeshListWithEdges);
      if (enableGLSLPhong)
        glslPhong->Enable();
    }
  }

  glStencilFunc(GL_ALWAYS, 0, ~(0u));

  // turns off stencil modifications
  glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

  glslPhong->Disable();

  if (renderExternalVertices)
  {
    glColor3f(0,1,0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for(int i=0; i<numExternalVertices; i++)
      glVertex3dv(&externalVertices[3*i]);
    glEnd();
  }

  // render current selected group
  if (currentSelectedObject >= 0 && currentSelectedGroup >= 0)
  {
    std::string groupNameStr = objMeshes[currentSelectedObject]->getGroupHandle(currentSelectedGroup)->getName();
    const char * temp = groupNameStr.c_str();
    char groupName[1024];
    strcpy(groupName, temp);
    // using color aqua
    glColor3f(0.0f, 1.0f, 1.0f);
    glLineWidth(5.0);
    objMeshRenders[currentSelectedObject]->renderGroupEdges(groupName);
  }

  if (enableVertexSelection)
  {
    glColor3f(1,0,0);
    glPointSize(6.5);
    glBegin(GL_POINTS);

    for(unsigned int objectIndex=0; objectIndex<objMeshes.size(); objectIndex++)
      if (selectedVerticesSet[objectIndex].size() != 0)
      {
        for(std::set<int>::iterator iter = selectedVerticesSet[objectIndex].begin(); iter != selectedVerticesSet[objectIndex].end(); iter++)
        {
          Vec3d v = objMeshes[objectIndex]->getPosition(*iter);
          glVertex3f(v[0], v[1], v[2]);
        }
      }
    glEnd();
  }

  if (renderLights)
  {
    if (lighting == NULL)
    {
      GLenum lightArray[8] = { 
        GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3,
        GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 };
  
      //GLfloat lP0[4] = { 0.2, 0.2, -0.2, 1.0 }; 
      GLfloat lP0[4] = { lightBoxMin[0], lightBoxMin[1], lightBoxMin[2], 1.0 };
      GLfloat lP1[4] = { lightBoxMax[0], lightBoxMin[1], lightBoxMin[2], 1.0 };
      GLfloat lP2[4] = { lightBoxMax[0], lightBoxMax[1], lightBoxMin[2], 1.0 };
      GLfloat lP3[4] = { lightBoxMin[0], lightBoxMax[1], lightBoxMin[2], 1.0 };
      GLfloat lP4[4] = { lightBoxMin[0], lightBoxMin[1], lightBoxMax[2], 1.0 };
      GLfloat lP5[4] = { lightBoxMax[0], lightBoxMin[1], lightBoxMax[2], 1.0 };
      GLfloat lP6[4] = { lightBoxMax[0], lightBoxMax[1], lightBoxMax[2], 1.0 };
      GLfloat lP7[4] = { lightBoxMin[0], lightBoxMax[1], lightBoxMax[2], 1.0 };
      GLfloat * lightPos[8] = { lP0, lP1, lP2, lP3, lP4, lP5, lP6, lP7 };
    
      glDisable(GL_LIGHTING);
      glColor3f(1,0,0);
      glPointSize(5.0);
      for(int light=0; light<8; light++)
      {
        if (glIsEnabled(lightArray[light]) == GL_TRUE)
        {
          glBegin(GL_POINTS);
            glVertex3f(lightPos[light][0], lightPos[light][1], lightPos[light][2]);
          glEnd();
          //printf("Light %d enabled at %G %G %G\n", light, lightPos[light][0], lightPos[light][1], lightPos[light][2]);
        }
      }
    }
    else
    {
      for(int light=0; light<8; light++)
      {
        if (lighting->IsLightEnabled(light))
        {
          const int LightPosArraySize = 4;
          float lightPos[LightPosArraySize];
          lighting->GetLightPosition(light, lightPos);
          glDisable(GL_LIGHTING);
          glColor3f(0,0,1);
          glPointSize(5.0);
          glBegin(GL_POINTS);
            glVertex3f(lightPos[0], lightPos[1], lightPos[2]);
          glEnd();
          //printf("Light %d enabled at %G %G %G\n", light, lightPos[0], lightPos[1], lightPos[2]);
        }
      }
    }
  }

  if (displayRuntimeInfo)
    DisplayRuntimeInformation();

  glutSwapBuffers();
}

Vec3d JetColorMap(double x)
{
  double a; // alpha
                                                                                                                                                             
  if(x < 0)
      return Vec3d(0,0,0);
  else if (x < 0.125) {
      a = x/0.125;
      return Vec3d(0, 0, 0.5+0.5*a);
  }
  else if (x < 0.375) {
      a = (x - 0.125)/0.25;
      return Vec3d(0, a, 1);
  }
  else if (x < 0.625) {
      a = (x - 0.375)/0.25;
      return Vec3d(a, 1, 1-a);
  }
  else if (x < 0.875) {
      a = (x - 0.625)/0.25;
      return Vec3d(1, 1-a, 0);
  }
  else if (x <= 1.0) {
      a = (x - 0.875)/0.125;
      return Vec3d(1-0.5*a, 0, 0);
  }
  else {
      return Vec3d(1,1,1);
  }
}


int main( int argc, char** argv )
{
  const int numFixedArguments = 2;
  if( argc < numFixedArguments )
  {
    std::cout << "Usage: " << argv[0] << " [obj file] [-f field filename] [-i] [-k] [-l<lighting file>] [-m] [-n] [-o<screenshot file>] [-s] [-v vertex file] [-w] [-0 file] [-1 file] [-y] [-z]" << std::endl;
    std::cout << "   -f: color vertices according to data from the given file (binary format file)" << std::endl;
    std::cout << "   -i: set light intensity"  << std::endl;
    std::cout << "   -l: set lighting filename"  << std::endl;
    std::cout << "   -L: first parameter is a list of obj files"  << std::endl;
    std::cout << "   -m: do not use lighting and materials" << std::endl;
    std::cout << "   -n: render using face normals" << std::endl;
    std::cout << "   -o: save a screen shot to file (support jpg, tif, ppm)" <<std::endl;
    std::cout << "   -s: use specular light component" << std::endl;
    std::cout << "   -t: only render in text to terminal" << std::endl;
    std::cout << "   -v: also render vertices from the given binary file (3 rows)" << std::endl;
    std::cout << "   -w: render edges" << std::endl;
    std::cout << "   -y: auto rotate model to make the longest bounding box edge along y-axis (works when only 1 object in the scene and no -z flag)" << std::endl;
    std::cout << "   -z: rotate z-up input model so that y-axis is up"<< std::endl;
    std::cout << "   -0: also render vertices from the given 0-indexed list of vertices" << std::endl;
    std::cout << "   -1: also render vertices from the given 1-indexed list of vertices" << std::endl;
    return 1;
  }

  objMeshname = argv[1];
  char fieldColorFilename[4096] = "__default";
  char lightIntensityString[4096] = "__default";
  char renderNormalLengthString[4096] = "__default";
  char renderExternalVerticesBinaryFilename[4096] = "__default";
  char renderExternalVerticesZeroIndexedFilename[4096] = "__default";
  char renderExternalVerticesOneIndexedFilename[4096] = "__default";
  renderMaterials = false;
  bool useFaceNormals = false;
  bool renderObjList = false;
  renderExternalVertices = false;
  bool renderToTerminal_ = false;

  opt_t opttable[] =
  {
    { (char*)"f", OPTSTR, fieldColorFilename },
    { (char*)"i", OPTSTR, &lightIntensityString },
    { (char*)"k", OPTINT, &numLights },
    { (char*)"l", OPTSTR, lightingConfigFilename },
    { (char*)"L", OPTBOOL, &renderObjList },
    { (char*)"m", OPTBOOL, &renderMaterials },
    { (char*)"n", OPTBOOL, &useFaceNormals },
    { (char*)"N", OPTSTR, &renderNormalLengthString },
    { (char*)"o", OPTSTR, &screenShotFilename },
    { (char*)"s", OPTBOOL, &showSpecularLight },
    { (char*)"t", OPTBOOL, &renderToTerminal_ },
    { (char*)"v", OPTSTR, &renderExternalVerticesBinaryFilename },
    { (char*)"w", OPTBOOL, &displayEdges },
    { (char*)"y", OPTBOOL, &autoRotateModelFlag },
    { (char*)"z", OPTBOOL, &rotateZUpToYUpFlag },
    { (char*)"0", OPTSTR, &renderExternalVerticesZeroIndexedFilename },
    { (char*)"1", OPTSTR, &renderExternalVerticesOneIndexedFilename },
    { NULL, 0, NULL }
  };
                                                                                                                                                             
  argv += numFixedArguments - 1;
  argc -= numFixedArguments - 1;
  int optup = getopts(argc,argv,opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %d: %s.\n",optup,argv[optup]);
    return 1;
  }

  renderMaterials = !renderMaterials;

  glutInit(&argc,argv);

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);

  //windowWidth = 1280;
  //windowHeight = 1024;

  windowWidth = 640;
  windowHeight = 640;

  if (strcmp(screenShotFilename, "__default") != 0)
  {
    windowWidth = SCREEN_SHOT_DEFAULT_WIDTH;
    windowHeight = SCREEN_SHOT_DEFAULT_HEIGHT;    
  }

  //windowWidth = 800;
  //windowHeight = 800;

  glutInitWindowSize (windowWidth,windowHeight);
  glutInitWindowPosition (0,0);
  winID = glutCreateWindow (objMeshname);

  //glutFullScreen();

  glutDisplayFunc(Display);

  glutReshapeFunc(Reshape);

  if (renderMaterials)
    GLUI_Master.set_glutIdleFunc( DoIdle );
  else
  {
    glutIdleFunc(DoIdle);
    enableGLSLPhong = 0;
    displayEdges = 1;
  }

  /* callback for mouse drags */
  glutMotionFunc(MouseDrag);

  /* callback for mouse movement without any buttons pressed */
  glutPassiveMotionFunc(MouseNoDrag);

  /* callback for mouse button changes */
  glutMouseFunc(MouseButtonActivity);

  glutKeyboardFunc(KeyboardFunc);

  glutSpecialFunc(SpecialKeysFunc);

  // parse input parameters
  if (strcmp(screenShotFilename,"__default") != 0)
  {
    // check if the file name is legal or not (in case, sometimes we input the obj filename by mistake)
    int strLength = (int)strlen(screenShotFilename);
    if (strLength <= 3)
    {
      printf("Error: screenshot filename %s is illegal.\n", screenShotFilename);
      exit(0);
    }
    char * filenameExt = &screenShotFilename[strlen(screenShotFilename)-4];
    bool filenameOKFlag = false;
    if (strcmp(filenameExt, ".jpg") == 0)
      filenameOKFlag = true;
    if (!filenameOKFlag && (strcmp(filenameExt, ".ppm") == 0))
      filenameOKFlag = true;
    if (!filenameOKFlag && (strcmp(filenameExt, ".tif") == 0))
      filenameOKFlag = true;
    
    if (!filenameOKFlag)
    {
      printf("Error: unsupported screenshot file type. (supported types are jpg, ppm, tif)\n");
      exit(0);
    }
  }

  if (strcmp(renderExternalVerticesBinaryFilename,"__default") != 0)
  {
    renderExternalVertices = true;
    int m1;
    ReadMatrixFromDisk_(renderExternalVerticesBinaryFilename,&m1,&numExternalVertices,&externalVertices);
    Assert_(m1,3,0);
  }

  if (strcmp(lightIntensityString,"__default") != 0)
    lightIntensity = strtod(lightIntensityString,NULL);

  if (strcmp(renderNormalLengthString,"__default") != 0)
    renderNormalLength = strtod(renderNormalLengthString,NULL);

  /* ############ generate string of input files ############ */ 
  std::vector<std::string> inputFiles;

  if (!renderObjList)
    inputFiles.push_back(objMeshname);
  else
  {
    FILE * fin;
    char mode[4] = "r";
    OpenFile_(objMeshname, &fin, mode);
    char s[4096];
    while(fgets(s,4096,fin) != NULL)
    {
      if (s[strlen(s)-1] == '\n')
        s[strlen(s)-1] = '\0';
      inputFiles.push_back(std::string(s));
    }
    fclose(fin);
  }

  /* ########### parse input files ################### */

  selectedVerticesSet = new std::set<int> [inputFiles.size()];
  selectedGroupsSet = new std::set<int> [inputFiles.size()];
  
  int numFaces=0;
  numTextures=0;
  for(unsigned int i=0; i<inputFiles.size(); i++)
  {
    printf("*** Loading %s (%d / %d)\n",inputFiles[i].c_str(),i+1,(int)inputFiles.size());
    ObjMesh * objMesh = new ObjMesh(inputFiles[i]);

    objMesh->buildFaceNormals();
    objMesh->buildVertexFaceNeighbors();
    objMesh->buildVertexNormals(85.0);

    numFaces += objMesh->getNumFaces();
    objMesh->getCubicBoundingBox(1.5,&lightBoxMin,&lightBoxMax); // 3.0
    printf("*** Bounding box size is %G\n",(lightBoxMax[0]-lightBoxMin[0])/3);

    if (useFaceNormals)
    {
      printf("*** Setting normals to face normals...\n");
      objMesh->setNormalsToFaceNormals();
    }
    objMeshes.insert(std::make_pair(i,objMesh));
    objMeshRenders.insert(std::make_pair(i, new ObjMeshRender(objMesh)));
    int textureMode = OBJMESHRENDER_GL_USEANISOTROPICFILTERING | OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_GL_MODULATE;
    //int textureMode = OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_GL_MODULATE;
    //int textureMode = OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_GL_REPLACE;
    //int textureMode = OBJMESHRENDER_GL_NOMIPMAP | OBJMESHRENDER_GL_MODULATE;
    int numFileTextures = objMeshRenders[i]->numTextures();
    numTextures += numFileTextures;
    if (numFileTextures > 0)
      objMeshRenders[i]->loadTextures(textureMode);
  }
  printf("*** Rendering %d faces.\n",numFaces);
  printf("*** Total num textures: %d.\n", numTextures);

  if (numTextures > 0)
  {
    printf("*** Disabling GLSL phong (due to textures).\n");
    enableGLSLPhong = 0;
  }

  if (rotateZUpToYUpFlag && autoRotateModelFlag)
  {
    printf("Warning: cannot support [-y] [-z] mode simultaneously. Auto-rotation [-y] functionality is disabled.\n");
    autoRotateModelFlag = false;
  }

  /* ################################################# */

  // ### load any vertices to be displayed ###
  
  if (strcmp(renderExternalVerticesZeroIndexedFilename,"__default") != 0)
  {
    // parse the 0-indexed list
    LoadList list;
    int * vertexIndices;
    list.load(renderExternalVerticesZeroIndexedFilename, &numExternalVertices, &vertexIndices);
    printf("Detected %d external vertices to render.\n", numExternalVertices);
    renderExternalVertices = true;
    externalVertices = (double*) malloc (sizeof(double) * 3 * numExternalVertices);
    for(int i=0; i<numExternalVertices; i++)
    {
      Vec3d pos = objMeshes[0]->getPosition(vertexIndices[i]);
      externalVertices[3*i+0] = pos[0];
      externalVertices[3*i+1] = pos[1];
      externalVertices[3*i+2] = pos[2];
      printf("0-indexed vtx %d: %G %G %G\n", vertexIndices[i], pos[0], pos[1], pos[2]);
    }
    free(vertexIndices);
  }

  if (strcmp(renderExternalVerticesOneIndexedFilename,"__default") != 0)
  {
    // parse the 1-indexed list
    LoadList list;
    int * vertexIndices;
    list.load(renderExternalVerticesOneIndexedFilename, &numExternalVertices, &vertexIndices);
    //list.printList(numExternalVertices, vertexIndices);
    printf("Detected %d external vertices to render.\n", numExternalVertices);
    renderExternalVertices = true;
    externalVertices = (double*) malloc (sizeof(double) * 3 * numExternalVertices);
    for(int i=0; i<numExternalVertices; i++)
    {
      Vec3d pos = objMeshes[0]->getPosition(vertexIndices[i]-1);
      externalVertices[3*i+0] = pos[0];
      externalVertices[3*i+1] = pos[1];
      externalVertices[3*i+2] = pos[2];
      printf("1-indexed vtx %d: %G %G %G\n", vertexIndices[i], pos[0], pos[1], pos[2]);
    }
    free(vertexIndices);
  }

  useFieldColors = false;
  if (strcmp(fieldColorFilename,"__default") != 0)
  {
    useFieldColors = true;

    if (renderObjList)
    {
      printf("Error: can't render field colors when rendering multiple input obj files.\n");
      return 1;
    }

    double * colors;
    int m1,n1;
    ReadMatrixFromDisk_(fieldColorFilename,&m1,&n1,&colors);
    Assert_(n1,1,0);
    Assert_(m1,(objMeshes[0])->getNumVertices(),1);

    // find max color
    double maxColor = -DBL_MAX;
    for(int i=0; i<m1; i++)
    {
      if (colors[i] > maxColor)
        maxColor = colors[i];
    }

    // find min color
    double minColor = DBL_MAX;
    for(int i=0; i<m1; i++)
    {
      if (colors[i] < minColor)
        minColor = colors[i];
    }

    double invMaxColor;
    if (maxColor <= 0)
      invMaxColor = 0;
    else
      invMaxColor = 1.0 / maxColor;

    for(int i=0; i<m1; i++)
    {
      double value = colors[i] * invMaxColor;
      fieldColors.push_back(JetColorMap(value));
    }

    printf("Max field value was: %G\n",maxColor);
    printf("Min field value was: %G\n",minColor);

    free(colors);
  }

  if (renderToTerminal_)
  {
    objMeshRenders[0]->outputOpenGLRenderCode();
    return 0;
  }

  Initialize();

  if (renderMaterials)
    InitGLUI();

  glutMainLoop();

  return(0);
}

