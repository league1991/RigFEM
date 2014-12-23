/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Reduced deformable dynamics" real-time driver application.           *
 * Uses model reduction to rapidly simulate deformable objects           *
 * undergoing large deformations.                                        *
 *                                                                       *
 * Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                            *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  Real-time reduced Saint-Venant Kirchhoff simulator.
  This is the demo code for the "StVK" SIGGRAPH 2005 paper:

  Jernej Barbic, Doug L. James: Real-Time Subspace Integration for 
  St.Venant-Kirchhoff Deformable Models, ACM Transactions on Graphics 24(3) 
  (SIGGRAPH 2005), p. 982-990, Los Angeles, CA, August 2005

  This code requires "GLUI", a LGPL-licensed GLUT-based UI library:
  http://glui.sourceforge.net/
  http://www.cs.unc.edu/~rademach/glui/
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cassert>
using namespace std;

#ifdef WIN32
  #include <windows.h>
#endif

#ifdef __APPLE__
  #include "TargetConditionals.h"
#endif

#include "initGraphics.h"
#include "sceneObjectReducedCPU.h"
#include "sceneObjectReducedGPU.h"
#include "performanceCounter.h"
#include "implicitNewmarkDense.h"
#include "implicitBackwardEulerDense.h"
#include "modalMatrix.h"
#include "reducedStVKForceModel.h"
#include "reducedLinearStVKForceModel.h"
#include "configFile.h"
#include "GL/glui.h"
#include "lighting.h"
#include "matrixIO.h"

char windowTitleBase[4096] = "Reduced StVK Demo";
void displayFunction(void);
void handleIdle(void);
int windowID;
int windowWidth = 800;
int windowHeight = 600;
double zNear=0.01;
double zFar=10.0;
double cameraRadius;
double focusPositionX, focusPositionY, focusPositionZ;
double cameraLongitude, cameraLattitude;
SphericalCamera * camera = NULL;
int g_iLeftMouseButton=0, g_iMiddleMouseButton=0, g_iRightMouseButton=0;
int g_vMousePos[2] = {0,0};
int shiftPressed=0;
int renderWireframe=1;
int renderAxes=0;
int renderDeformableObject=1;
int enableNormalCorrection = 0;
int computeDynamicNormals = 0;
int lockScene=0;
int staticSolver = 0;
int renderOnGPU = 1;
int displayContactInfo = 0;
int syncTimeStepWithGraphics=1;
float timeStep = 1.0 / 30;
float newmarkBeta = 0.25;
float newmarkGamma = 0.5;
int use1DNewmarkParameterFamily = 1;
int substepsPerTimeStep = 1;
double fps = 30.0;
double cpuLoad = 0;
int enableTextures = 0;
int forceModel = 0;
double * uClosestVertexBuffer = NULL;

int graphicFrame = 0;
int pulledVertex = -1;
int dragStartX, dragStartY;
int explosionFlag = 0;

vector<string> configFilenames;

PerformanceCounter titleBarCounter;
PerformanceCounter explosionCounter;
Lighting * lighting = NULL;
SceneObjectReduced * deformableObjectRenderingMeshReduced = NULL;
SceneObjectReducedGPU * deformableObjectRenderingMeshGPU = NULL;
SceneObjectReducedCPU * deformableObjectRenderingMeshCPU = NULL;
SceneObject * extraSceneGeometry = NULL;
ModalMatrix * renderingModalMatrix = NULL;
ImplicitNewmarkDense * implicitNewmarkDense = NULL;
StVKReducedInternalForces * stVKReducedInternalForces = NULL;
StVKReducedStiffnessMatrix * stVKReducedStiffnessMatrix = NULL;
ReducedForceModel * reducedForceModel = NULL;
ReducedStVKForceModel * reducedStVKForceModel;
ReducedLinearStVKForceModel * reducedLinearStVKForceModel;
double lowestFrequency;

int plasticDeformationsEnabled = 0;
float plasticThreshold = 1E9;

void stopDeformations_buttonCallBack(int code);

PerformanceCounter cpuLoadCounter;

int nRendering;
int r;
double * q = NULL;
double * fq = NULL;
double * fqBase = NULL;

// options for the config file
char deformableObjectFilename[4096];
char extraSceneGeometryFilename[4096];
char modesFilename[4096];
char cubicPolynomialFilename[4096];
char lightingConfigFilename[4096];
float dampingMassCoef;
float dampingStiffnessCoef;
char backgroundColorString[4096] = "255 255 255";

GLUI * glui;
GLUI_Spinner * timeStep_spinner;
GLUI_Checkbox * renderOnGPU_checkbox;
float deformableObjectCompliance = 1.0;
float frequencyScaling = 1.0;
int sceneID=0;

void callAllUICallBacks();
void cleanupScene();
void Sync_GLUI();

int developerMode = 0;
string configFilename;
string configFilesDir;

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

// called periodically by GLUT (the main graphics function)
void displayFunction(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  // Setup model transformations.
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();

  camera->Look();

  // if GPU rendering, must set lighting 
  // else set lights directly via OpenGL (all done inside the following function)
  deformableObjectRenderingMeshReduced->SetLighting(lighting);

  // when GPU-rendering, set the same lighting for the CPU-drawn objects
  if (renderOnGPU)
    lighting->LightScene();

  glEnable(GL_LIGHTING);

  // the stencil buffer is used to determine the 3D vertex closest to a mouse click
  glStencilFunc(GL_ALWAYS, 1, ~(0u));

  // render the deformable object
  if (renderDeformableObject)
  {
    deformableObjectRenderingMeshReduced->Render();
  }

  glStencilFunc(GL_ALWAYS, 0, ~(0u));

  // render any extra scene geometry
  if (extraSceneGeometry != NULL)
    extraSceneGeometry->Render();

  glStencilFunc(GL_ALWAYS, 1, ~(0u));

  glDisable(GL_LIGHTING);

  // render the deformable object wireframe
  glColor3f(0,0,0);
  if (renderWireframe)
  {
    deformableObjectRenderingMeshReduced->RenderEdges();
  }

  glStencilFunc(GL_ALWAYS, 0, ~(0u));

  glColor3f(0,0,0);

  if (renderAxes)
    drawAxes(1.0);

  // render the currently pulled vertex
  if (pulledVertex >= 0)
  {
    double pulledVertexPos[3];
    deformableObjectRenderingMeshReduced->GetSingleVertexPosition(pulledVertex,
      &pulledVertexPos[0], &pulledVertexPos[1], &pulledVertexPos[2]);
    glColor3f(0,1,0);

    glEnable(GL_POLYGON_OFFSET_POINT);
    glPolygonOffset(-1.0,-1.0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
      glVertex3f(pulledVertexPos[0], pulledVertexPos[1], pulledVertexPos[2]);
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
  }

  // ==== bitmap routines below here
  glMatrixMode(GL_MODELVIEW); 
  glPushMatrix(); 
  glLoadIdentity(); 
  glMatrixMode(GL_PROJECTION); 
  glPushMatrix(); 
  glLoadIdentity(); 

  // print info in case of integrator blow-up
  char s[4096];
  glColor3f(1,0,0);
  if (explosionFlag)
  {
    sprintf(s,"The integrator went unstable.");
    double x1 = 10;
    double y1 = 25;
    double X1 = -1 + 2.0 * x1 / windowWidth;
    double Y1 = -1 + 2.0 * y1 / windowHeight;
    print_bitmap_string(X1,Y1,-1,GLUT_BITMAP_9_BY_15 ,s);

    sprintf(s,"Reduce the timestep, or increase the number of substeps per timestep.");
    x1 = 10;
    y1 = 10;
    X1 = -1 + 2.0 * x1 / windowWidth;
    Y1 = -1 + 2.0 * y1 / windowHeight;
    print_bitmap_string(X1,Y1,-1,GLUT_BITMAP_9_BY_15 ,s);
  }

  glPopMatrix(); 
  glMatrixMode(GL_MODELVIEW); 
  glPopMatrix(); 

  glutSwapBuffers();
}
                                
// the "idle" routine; called periodically by GLUT 
void idleFunction(void)
{
  cpuLoadCounter.StartCounter();

  glutSetWindow(windowID);
  
  if (!lockScene)
  {
    // determine force in case user is pulling on a vertex
    if (g_iLeftMouseButton) 
    {
      if (pulledVertex != -1)
      {
        double forceX = (g_vMousePos[0] - dragStartX);
        double forceY = -(g_vMousePos[1] - dragStartY);

        double externalForce[3];

        camera->CameraVector2WorldVector_OrientationOnly3D(
          forceX, forceY, 0, externalForce);

        renderingModalMatrix->ProjectSingleVertex(pulledVertex,
          externalForce[0], externalForce[1], externalForce[2], fq);

        for(int i=0; i<r; i++)
          fq[i] = fqBase[i] + deformableObjectCompliance * fq[i];
      }
    }
    else
    {
      memcpy(fq,fqBase,sizeof(double) * r);
    }

    // set the reduced external forces
    implicitNewmarkDense->SetExternalForces(fq);

    // integrate the dynamics via implicit Newmark
    for(int i=0; i<substepsPerTimeStep; i++)
    {
      int code = implicitNewmarkDense->DoTimestep();
      if (code != 0)
      {
        printf("The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep.\n");
        implicitNewmarkDense->ResetToRest();
        for(int i=0; i<r; i++)
        {
          fqBase[i] = 0;
          fq[i] = 0;
        }
        implicitNewmarkDense->SetExternalForces(fq);
        explosionFlag = 1;
        explosionCounter.StartCounter();
        break;
      }

      /*
        printf("q =\n");
        double * q = implicitNewmarkDense->Getq();
        for(int i=0; i<r; i++)
          printf("%G ", q[i]);
        printf("\n");
      */
    }

    memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);
  }

  if (explosionFlag)
  {
    explosionCounter.StopCounter();
    if (explosionCounter.GetElapsedTime() > 4.0) // the message will appear on screen for 4 seconds
      explosionFlag = 0;
  }

  // compute u=Uq
  deformableObjectRenderingMeshReduced->Setq(q);
  deformableObjectRenderingMeshReduced->Compute_uUq();

  graphicFrame++;
  
  // update title bar information at 4 Hz
  titleBarCounter.StopCounter();
  double elapsedTime = titleBarCounter.GetElapsedTime();
  if (elapsedTime >= 1.0 / 4)
  {
    titleBarCounter.StartCounter();
    fps = graphicFrame / elapsedTime;

    // update menu bar
    char windowTitle[4096];
    sprintf(windowTitle,"%s | Num modes = %d | %.1f Hz | Deformation CPU Load: %d%%", windowTitleBase, 
      implicitNewmarkDense->GetNumDOFs() , fps, (int)(100 * cpuLoad + 0.5) );
    glutSetWindowTitle(windowTitle);
    graphicFrame = 0;

    if (syncTimeStepWithGraphics)
    {
      timeStep = 1.0 / fps;
      implicitNewmarkDense->SetTimestep(timeStep / substepsPerTimeStep);
      Sync_GLUI();
    }
  }

  cpuLoadCounter.StopCounter();
  double cpuTimePerGraphicsFrame = cpuLoadCounter.GetElapsedTime();
  cpuLoad = cpuTimePerGraphicsFrame * fps; 

  glutPostRedisplay();
}

void keyboardFunction (unsigned char key, int x, int y)
{
  double cameraX,cameraY,cameraZ;

  switch (key)
  {
    case 27:
      exit(0);
    break;

    case 'w':
      renderWireframe = !renderWireframe;
    break;

    case 'i':
      camera->GetAbsWorldPosition(cameraX,cameraY,cameraZ);
      printf("Camera is positioned at: %G %G %G\n",cameraX,cameraY,cameraZ);
      printf("Camera radius is: %G \n",camera->GetRadius());
      printf("Camera Phi is: %G \n",180.0/M_PI*camera->GetPhi());
      printf("Camera Theta is: %G \n",180.0/M_PI*camera->GetTheta());
    break;

    case '\\':
      camera->Reset();
    break;

    case 'a':
      renderAxes = !renderAxes;
    break;

    // in this mode, can move the camera while the object's deformations are frozen
    case 'l':
      lockScene = !lockScene;
      if (lockScene)
      {
        camera->PushPosition();
      }
      else
      {
        camera->PopPosition();
      }
    break;

    case 'e':
      if (developerMode)
        renderDeformableObject = !renderDeformableObject;
    break;

    // print out reduced coordinates
    case 'q':
      for(int i=0; i<r; i++)
        printf("%G ", q[i]);
      printf("\n");
      printf("Also writing reduced coordinates to out.q .\n");
      WriteMatrixToDisk("out.q", r,1, q);
    break;

    // make the current configuration a new rest configuration
    case 'z':
      for(int i=0; i<r; i++)
        fqBase[i] = fq[i];
      dragStartX = g_vMousePos[0];
      dragStartY = g_vMousePos[1];
      break;

    // reset the rest configuration to zero deformation
    case 'Z':
      for(int i=0; i<r; i++)
        fqBase[i] = 0;
      dragStartX = g_vMousePos[0];
      dragStartY = g_vMousePos[1];
      break;

    case '0':
      stopDeformations_buttonCallBack(0);
      break;
  }
}

void reshape(int x, int y)
{
  glViewport(0,0,x,y);

  windowWidth = x;
  windowHeight = y;

  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 

  gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, zNear, zFar);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void exitHandler()
{
  printf("Executing the exit handler...\n");
  printf("De-allocating data structures...\n");
  cleanupScene();
}

void mouseMotionFunction(int x, int y)
{
  int mouseDeltaX = x-g_vMousePos[0];
  int mouseDeltaY = y-g_vMousePos[1];

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;

  if (g_iLeftMouseButton) // left mouse button only responds to clicks (not to drag)
  {
  }

  if (g_iRightMouseButton) // right mouse button handles camera rotations
  {
    const double factor = 0.2;
    camera->MoveRight(factor * mouseDeltaX);
    camera->MoveUp(factor * mouseDeltaY);
  }

  if (g_iMiddleMouseButton) // handle zoom in/out
  {
    const double factor = 0.1;
    camera->ZoomIn(cameraRadius * factor * mouseDeltaY);
  }
}

void mouseButtonActivityFunction(int button, int state, int x, int y)
{
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      shiftPressed = g_iLeftMouseButton && (glutGetModifiers() == GLUT_ACTIVE_SHIFT);

      if (g_iLeftMouseButton) 
      {
        // check if user clicked on the object, if yes, determine closest 3D vertex to the click location

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

        if (stencilValue == 1)
        {
          dragStartX = x;
          dragStartY = y;
          Vec3d queryPos(worldX, worldY, worldZ);
          if (renderOnGPU)
            pulledVertex = deformableObjectRenderingMeshReduced->GetClosestVertex(
              queryPos, NULL, uClosestVertexBuffer);
          else
          {
            pulledVertex = deformableObjectRenderingMeshReduced->GetClosestVertex(queryPos, NULL, NULL);
          }
          printf("Clicked on vertex: %d\n", pulledVertex);
        }
      }

      if (!g_iLeftMouseButton)
      {
        pulledVertex = -1;
      }

    break;

    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      shiftPressed = g_iLeftMouseButton && (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
    break;

    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      shiftPressed = g_iLeftMouseButton && (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
    break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void initScene()
{
  // init lighting
  try
  {
    lighting = new Lighting(lightingConfigFilename);
  }
  catch(int exceptionCode)
  {
    printf("Error (%d) reading lighting information from %s .\n", exceptionCode, lightingConfigFilename);
    exit(1);
  }

  // init camera
  delete(camera);
  double virtualToPhysicalPositionFactor = 1.0;
  initCamera(cameraRadius, cameraLongitude, cameraLattitude, focusPositionX, focusPositionY, focusPositionZ, 1.0 / virtualToPhysicalPositionFactor, &zNear, &zFar, &camera);

  // load the rendering modes of the deformable object
  int nRendering;
  float * URenderingFloat;
  ReadMatrixFromDisk_(modesFilename,&nRendering,&r,&URenderingFloat);
  nRendering /= 3;
  double * URendering = (double*) malloc (sizeof(double) * 3 * nRendering * r);
  for(int i=0; i < 3 * nRendering * r; i++)
    URendering[i] = URenderingFloat[i];
  free(URenderingFloat);
  renderingModalMatrix = new ModalMatrix(nRendering,r,URendering);
  free(URendering); // ModalMatrix made an internal copy

  // init room for reduced coordinates and reduced forces
  q = (double*) calloc (r, sizeof(double));
  fq = (double*) calloc (r, sizeof(double));
  fqBase = (double*) calloc (r, sizeof(double));

  // initialize the GPU rendering class for the deformable object
  try
  {
    deformableObjectRenderingMeshGPU = new SceneObjectReducedGPU(deformableObjectFilename, renderingModalMatrix); // uses GPU to compute u=Uq
  }
  catch(int exceptionCode)
  {
    printf("Warning: unable to initialize GPU rendering for deformations (code: %d). Using CPU instead.\n", exceptionCode);
    deformableObjectRenderingMeshGPU = NULL;
    renderOnGPU = 0;
  }

  // initialize the CPU rendering class for the deformable object
  deformableObjectRenderingMeshCPU = new SceneObjectReducedCPU(deformableObjectFilename, renderingModalMatrix); // uses CPU to compute u=Uq

  // prepare textures (if necessary)
  if (enableTextures)
  {
    if (deformableObjectRenderingMeshGPU != NULL)
    {
      deformableObjectRenderingMeshGPU->SetUpTextures(SceneObject::MODULATE, SceneObject::USEMIPMAP);
      deformableObjectRenderingMeshGPU->EnableTextures();
    }
    deformableObjectRenderingMeshCPU->SetUpTextures(SceneObject::MODULATE, SceneObject::USEMIPMAP);
    deformableObjectRenderingMeshCPU->EnableTextures();
  }

  // create a buffer used when searching for the 3D vertex closest to the click location
  if (deformableObjectRenderingMeshGPU != NULL)
    uClosestVertexBuffer = (double*) 
      malloc (sizeof(double) * 3 * deformableObjectRenderingMeshGPU->Getn());

  if (renderOnGPU)
    deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshGPU;
  else
    deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshCPU;

  deformableObjectRenderingMeshReduced->ResetDeformationToRest();

  // make reduced mass matrix (=identity)
  double * massMatrix = (double*) calloc (r*r,sizeof(double));
  for(int i=0; i<r; i++)
    massMatrix[ELT(r,i,i)] = 1.0;

  // load internal force cubic polynomial
  #if defined(__APPLE__) && (TARGET_CPU_PPC)
    stVKReducedInternalForces = new StVKReducedInternalForces(cubicPolynomialFilename, -1, 1); // version that can load a little endian binary file to a big endian machine
  #else
    stVKReducedInternalForces = new StVKReducedInternalForces(cubicPolynomialFilename); // normal version
  #endif

  // create stiffness matrix polynomials
  stVKReducedStiffnessMatrix = new StVKReducedStiffnessMatrix(stVKReducedInternalForces);

  // create the "internal force models" 
  reducedStVKForceModel = new ReducedStVKForceModel(stVKReducedInternalForces, stVKReducedStiffnessMatrix);
  reducedLinearStVKForceModel = new ReducedLinearStVKForceModel(stVKReducedStiffnessMatrix);
  reducedForceModel = reducedStVKForceModel;

  // init the implicit Newmark
  implicitNewmarkDense = new ImplicitNewmarkDense(r, timeStep, massMatrix, reducedForceModel, ImplicitNewmarkDense::positiveDefiniteMatrixSolver, dampingMassCoef, dampingStiffnessCoef);
  implicitNewmarkDense->SetTimestep(timeStep / substepsPerTimeStep);
  implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
  implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);
  free(massMatrix);

  // load any external geometry file (e.g. some static scene for decoration; usually there will be none)
  if (strcmp(extraSceneGeometryFilename,"__none") != 0)
  {
    extraSceneGeometry = new SceneObject(extraSceneGeometryFilename);
    extraSceneGeometry->BuildNormals(85.0);
  }
  else
    extraSceneGeometry = NULL;

  // compute lowest frequency of the system (smallest eigenvalue of K)
  double * K = (double*) malloc (sizeof(double) * r * r);
  double * zero = (double*) calloc (r, sizeof(double));
  stVKReducedStiffnessMatrix->Evaluate(zero, K);

/*
  // find smallest eigenvalue of K
  Matrix<double> KM(r, r, K, false, false);
  Matrix<double> EigenVectors(r,r);
  Matrix<double> Lambda(r,1);
  KM.SymmetricEigenDecomposition(EigenVectors, Lambda);
  lowestFrequency = Lambda(0,0);
  if (lowestFrequency < 0)
  {
    printf("Warning: negative eigenvalue encountered.\n");
    lowestFrequency = 1.0;
  }
  else
  {
    lowestFrequency = sqrt(lowestFrequency) / (2 * M_PI);
    printf("System lowest frequency is: %G\n", lowestFrequency);
  }
  free(zero);
  free(K);
*/

  // set background color
  int colorR, colorG, colorB;
  sscanf(backgroundColorString, "%d %d %d", &colorR, &colorG, &colorB);
  glClearColor(1.0 * colorR / 255, 1.0 * colorG / 255, 1.0 * colorB / 255, 0.0);

  Sync_GLUI();
  callAllUICallBacks();

}

void cleanupScene()
{
}

void initConfigurations()
{
  ConfigFile configFile;

  // specify the entries of the config file
  configFile.addOptionOptional("windowWidth",&windowWidth,800);
  configFile.addOptionOptional("windowHeight",&windowHeight,800);

  configFile.addOption("deformableObjectFilename",deformableObjectFilename);
  configFile.addOptionOptional("modesFilename",modesFilename,"__none");
  configFile.addOptionOptional("cubicPolynomialFilename",cubicPolynomialFilename,"__none");

  configFile.addOption("dampingMassCoef",&dampingMassCoef);
  configFile.addOption("dampingStiffnessCoef",&dampingStiffnessCoef);

  configFile.addOptionOptional("plasticThreshold", &plasticThreshold, plasticThreshold);

  configFile.addOption("deformableObjectCompliance",&deformableObjectCompliance);
  configFile.addOption("frequencyScaling",&frequencyScaling);

  configFile.addOptionOptional("cameraRadius",&cameraRadius,17.5);
  configFile.addOptionOptional("focusPositionX",&focusPositionX,0.0);
  configFile.addOptionOptional("focusPositionY",&focusPositionY,0.0);
  configFile.addOptionOptional("focusPositionZ",&focusPositionZ,0.0);
  configFile.addOptionOptional("cameraLongitude",&cameraLongitude,-60.0);
  configFile.addOptionOptional("cameraLattitude",&cameraLattitude,20.0);

  configFile.addOptionOptional("renderWireframe",&renderWireframe,1);

  configFile.addOptionOptional("extraSceneGeometry",extraSceneGeometryFilename,"__none");

  configFile.addOptionOptional("enableTextures",&enableTextures,enableTextures);

  configFile.addOptionOptional("backgroundColor",backgroundColorString, backgroundColorString);

  string lightingConfigFilenameDefault = configFilesDir + "default.lighting";
  configFile.addOptionOptional("lightingConfigFilename",lightingConfigFilename,
    (char*) lightingConfigFilenameDefault.c_str());

  configFile.addOptionOptional("substepsPerTimeStep", 
    &substepsPerTimeStep, substepsPerTimeStep);

  configFile.addOptionOptional("renderOnGPU", &renderOnGPU, 1);

  // parse the configuration file
  if (configFile.parseOptions((char*)configFilename.c_str()) != 0)
  {
    printf("Error: unable to load the configuration file.\n");
    exit(1);
  }
  // the config variables have now been loaded with their specified values

  // informatively print the variables (with assigned values) that were just parsed
  configFile.printOptions();
}

void Sync_GLUI()
{
  glui->sync_live();
}

void displayContactInfoCallBack(int code)
{
}

void sceneListBoxCallback(int code)
{
  cleanupScene();
  initConfigurations();
  initScene();
}

void deformableObjectCompliance_spinnerCallBack(int code)
{
  if (deformableObjectCompliance < 0)
    deformableObjectCompliance = 0;

  glui->sync_live();
}

void timeStep_spinnerCallBack(int code)
{
  if (timeStep < 0)
    timeStep = 0;

  implicitNewmarkDense->SetTimestep(timeStep / substepsPerTimeStep);

  glui->sync_live();
}

void syncTimeStepWithGraphics_checkboxCallBack(int code)
{
  if (syncTimeStepWithGraphics)
    timeStep_spinner->disable();
  else
    timeStep_spinner->enable();
}

void frequencyScaling_spinnerCallBack(int code)
{
  if (frequencyScaling < 0)
    frequencyScaling = 0;

  glui->sync_live();

  implicitNewmarkDense->SetInternalForceScalingFactor(frequencyScaling * frequencyScaling);
}

void newmarkBeta_spinnerCallBack(int code)
{
  if (newmarkBeta < 0)
    newmarkBeta = 0;

  if (newmarkBeta > 0.5)
    newmarkBeta = 0.5;

  if (use1DNewmarkParameterFamily)
  {
    if (newmarkBeta > 0.25)
      newmarkGamma = sqrt(4.0 * newmarkBeta) - 0.5;
    else
      newmarkGamma = 0.5;
  }

  implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
  implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);

  glui->sync_live();
}

void newmarkGamma_spinnerCallBack(int code)
{
  if (newmarkGamma < 0.5)
    newmarkGamma = 0.5;

  if (newmarkGamma > 1.0)
    newmarkGamma = 1.0;

  if (use1DNewmarkParameterFamily)
    newmarkBeta = (newmarkGamma + 0.5) * (newmarkGamma + 0.5) / 4.0;

  implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
  implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);

  glui->sync_live();
}

void newmark_checkboxuse1DNewmarkParameterFamilyCallBack(int code)
{
  if (use1DNewmarkParameterFamily)
  {
    newmarkBeta = (newmarkGamma + 0.5) * (newmarkGamma + 0.5) / 4.0;
    implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
    implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);
  }

  glui->sync_live();
}


void rayleighMass_spinnerCallBack(int code)
{
  if (dampingMassCoef < 0)
    dampingMassCoef = 0;

  implicitNewmarkDense->SetDampingMassCoef(dampingMassCoef);

  glui->sync_live();
}

void rayleighStiffness_spinnerCallBack(int code)
{
  if (dampingStiffnessCoef < 0)
    dampingStiffnessCoef = 0;

  implicitNewmarkDense->SetDampingStiffnessCoef(dampingStiffnessCoef);

  glui->sync_live();
}

void timeStepSubdivisions_spinnerCallBack(int code)
{
  if (substepsPerTimeStep < 1)
    substepsPerTimeStep = 1;

  implicitNewmarkDense->SetTimestep(timeStep / substepsPerTimeStep); 

  glui->sync_live();
}

void plasticDeformationsEnabled_checkboxCallBack(int code)
{
  implicitNewmarkDense->SetPlasticThreshold(plasticThreshold);
  implicitNewmarkDense->UsePlasticDeformations(plasticDeformationsEnabled);

  glui->sync_live();
}

void renderOnGPU_checkBoxCallBack(int code)
{
  if (renderOnGPU && (deformableObjectRenderingMeshGPU != NULL))
  {
    deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshGPU;
  }
  else
  {
    deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshCPU;
  }
}

void forceModel_checkBoxCallBack(int code)
{
  if (reducedForceModel == reducedStVKForceModel)
    reducedForceModel = reducedLinearStVKForceModel;
  else
    reducedForceModel = reducedStVKForceModel;
  implicitNewmarkDense->SetReducedForceModel(reducedForceModel);
}

void stopDeformations_buttonCallBack(int code)
{
  implicitNewmarkDense->ResetToRest();
  plasticDeformationsEnabled_checkboxCallBack(0);
}

void staticSolver_checkboxCallBack(int code)
{
  implicitNewmarkDense->UseStaticSolver(staticSolver);
}

void exit_buttonCallBack(int code)
{
  exit(0);
}

// calls all GLUI callbacks, except the listBox callbacks
void callAllUICallBacks()
{
  deformableObjectCompliance_spinnerCallBack(0);
  frequencyScaling_spinnerCallBack(0);
  timeStep_spinnerCallBack(0);
  syncTimeStepWithGraphics_checkboxCallBack(0);
  rayleighMass_spinnerCallBack(0);
  rayleighStiffness_spinnerCallBack(0);
  timeStepSubdivisions_spinnerCallBack(0);
  renderOnGPU_checkBoxCallBack(0);
  newmarkBeta_spinnerCallBack(0);
  newmarkGamma_spinnerCallBack(0);
  newmark_checkboxuse1DNewmarkParameterFamilyCallBack(0);
}

void initGLUI()
{
  // generate the UI

  glui = GLUI_Master.create_glui("Controls", 0, windowWidth + 10, 0);

  glui->add_spinner("Deformable object compliance:", GLUI_SPINNER_FLOAT, &deformableObjectCompliance, 0, deformableObjectCompliance_spinnerCallBack );

  glui->add_spinner("Frequency scaling:", GLUI_SPINNER_FLOAT, &frequencyScaling, 0, frequencyScaling_spinnerCallBack);

  // ******** newmark beta, gamma ********
  
  GLUI_Panel * newmark_panel =
    glui->add_panel("Newmark integrator parameters", GLUI_PANEL_EMBOSSED);
  newmark_panel->set_alignment(GLUI_ALIGN_LEFT);
  
  glui->add_checkbox_to_panel(newmark_panel, "Link Beta and Gamma", 
    &use1DNewmarkParameterFamily, 0, newmark_checkboxuse1DNewmarkParameterFamilyCallBack);

  GLUI_Spinner * newmarkBeta_spinner = 
    glui->add_spinner_to_panel(newmark_panel,"Beta", GLUI_SPINNER_FLOAT,
      &newmarkBeta, 0, newmarkBeta_spinnerCallBack);
  newmarkBeta_spinner->set_speed(0.1);

  GLUI_Spinner * newmarkGamma_spinner = 
    glui->add_spinner_to_panel(newmark_panel,"Gamma", GLUI_SPINNER_FLOAT,
      &newmarkGamma, 0, newmarkGamma_spinnerCallBack);
  newmarkGamma_spinner->set_speed(0.1);

    glui->add_checkbox_to_panel(newmark_panel,"Static solver only", &staticSolver, 0, staticSolver_checkboxCallBack);

  glui->add_checkbox_to_panel(newmark_panel, "Linear reduced model", &forceModel, 0, forceModel_checkBoxCallBack);

  // ******** damping ********

  GLUI_Panel * damping_panel =
    glui->add_panel("Tangential Rayleigh Damping", GLUI_PANEL_EMBOSSED);
  damping_panel->set_alignment(GLUI_ALIGN_LEFT);

  glui->add_spinner_to_panel(damping_panel,"Mass-proportional", GLUI_SPINNER_FLOAT,
      &dampingMassCoef, 0, rayleighMass_spinnerCallBack);

  glui->add_spinner_to_panel(damping_panel,"Stiffness-proportional", GLUI_SPINNER_FLOAT,
      &dampingStiffnessCoef, 0, rayleighStiffness_spinnerCallBack);

  glui->add_button("Stop deformations", 0, stopDeformations_buttonCallBack);

  // ******** timestep control ********

  GLUI_Panel * timeStep_panel =
    glui->add_panel("Timestep control", GLUI_PANEL_EMBOSSED);
  timeStep_panel->set_alignment(GLUI_ALIGN_LEFT);

  glui->add_checkbox_to_panel(timeStep_panel, "Sync with graphics", 
    &syncTimeStepWithGraphics, 0, syncTimeStepWithGraphics_checkboxCallBack);

  timeStep_spinner = 
    glui->add_spinner_to_panel(timeStep_panel,"Timestep [sec]", GLUI_SPINNER_FLOAT,
      &timeStep, 0, timeStep_spinnerCallBack);
  timeStep_spinner->set_alignment(GLUI_ALIGN_LEFT);

  if (syncTimeStepWithGraphics)
    timeStep_spinner->disable();

  glui->add_spinner_to_panel(timeStep_panel,"Substeps per timestep", GLUI_SPINNER_INT, &substepsPerTimeStep, 0, timeStepSubdivisions_spinnerCallBack);

  // ******* plastic deformations ********

  GLUI_Panel * plasticDeformationsPanel = glui->add_panel("Plastic deformations", GLUI_PANEL_EMBOSSED);
  plasticDeformationsPanel->set_alignment(GLUI_ALIGN_LEFT);
  glui->add_checkbox_to_panel(plasticDeformationsPanel, "Enable", &plasticDeformationsEnabled, 0, plasticDeformationsEnabled_checkboxCallBack);
  glui->add_column_to_panel(plasticDeformationsPanel, 0);
  glui->add_spinner_to_panel(plasticDeformationsPanel, "Threshold:", GLUI_SPINNER_FLOAT, &plasticThreshold, 0, plasticDeformationsEnabled_checkboxCallBack);

  // *** u = Uq ***

  renderOnGPU_checkbox = glui->add_checkbox("Compute u=Uq on GPU", &renderOnGPU, 0, renderOnGPU_checkBoxCallBack);

  glui->add_separator();

  GLUI_Panel * instructions_panel = glui->add_panel("Mouse buttons:", GLUI_PANEL_EMBOSSED);
  instructions_panel->set_alignment(GLUI_ALIGN_LEFT);
  glui->add_statictext_to_panel(instructions_panel, "Left + drag: apply force");
  glui->add_statictext_to_panel(instructions_panel, "Middle + drag: zoom in/out");
  glui->add_statictext_to_panel(instructions_panel, "Right + drag: rotate camera");
 
  glui->add_separator();

  glui->add_statictext("Jernej Barbic and Doug James");
  glui->add_statictext("Carnegie Mellon University, Cornell, 2007");

  glui->add_separator();

  glui->add_button("Exit program", 0, exit_buttonCallBack);

  Sync_GLUI();

  glui->set_main_gfx_window( windowID );
}

int main(int argc, char* argv[])
{
  if ( argc < 2 ) 
  {
    printf("Reduced StVK simulation of a large-deformation object.\n");
    printf("Usage: %s [config file]\n", argv[0]);
    return 1;
  }

  printf("Starting application.\n");
  configFilesDir = string("configFiles\\");
  configFilename = string(argv[1]);
  printf("Loading scene configuration from %s.\n", configFilename.c_str());

  initConfigurations();

  atexit(exitHandler);

  initGLUT(argc, argv, windowTitleBase , windowWidth, windowHeight, &windowID);
  initGraphics(windowWidth, windowHeight);

  initGLUI();

  initScene();

  // disable the GPU UI switch if no GPU support for vertex texture fetches
  if (deformableObjectRenderingMeshGPU == NULL)
    renderOnGPU_checkbox->disable();

  glutMainLoop(); 

  return 0;
}

