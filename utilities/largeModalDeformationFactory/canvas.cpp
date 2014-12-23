/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
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

#ifdef WIN32
  #include <windows.h>
#endif

#include "openGL-headers.h"
#include "objMeshRender.h"
#include "canvas.h"
#include "renderVolumetricMesh.h"

#ifndef M_PI
  #define M_PI 3.1415926535897932384
#endif

BEGIN_EVENT_TABLE(MyGLCanvas, wxGLCanvas)
  EVT_SIZE(MyGLCanvas::OnSize)
  EVT_PAINT(MyGLCanvas::OnPaint)
  EVT_ERASE_BACKGROUND(MyGLCanvas::OnEraseBackground)
  EVT_MOUSE_EVENTS(MyGLCanvas::OnMouse)
  EVT_KEY_DOWN(MyGLCanvas::OnKey)
  EVT_IDLE(MyGLCanvas::OnIdle)
END_EVENT_TABLE()

MyGLCanvas::MyGLCanvas(PrecomputationState * precomputationState, UIState * uiState, MyFrame * parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name, int* attribList, const wxPalette& palette) :
  wxGLCanvas(parent, id, attribList, pos, size, style | wxFULL_REPAINT_ON_RESIZE, name, palette) 
{
  m_context = new wxGLContext(this);
  this->parent = parent;
  this->precomputationState = precomputationState;
  this->uiState = uiState;

  cameraRadius = 5.0;

  double centerV[3] = {0, 0, 0};
  double upVector[3] = {0,1,0};

  // init camera for 3d viewing
  camera = new SphericalCamera(cameraRadius, 1.0 * 110 / 360 * (2*M_PI), 
    1.0 * 30 / 360 * (2*M_PI),
    centerV, upVector,  0.02, 1);

  renderingMeshDisplayListAvailable = false;
  showRenderingMesh = false;

  simulationMeshDisplayListAvailable = false;
  showSimulationMesh = false;

  showFixedVertices = false;

  showLinearModes = false;
  uLinear = NULL;

  showNonLinearModes = false;
  uNonLinear = NULL;

  renderedLinearMode = 0;
  linearRenderingMagnitude = 100.0;

  renderedNonLinearMode = 0;
  nonLinearRenderingMagnitude = 100.0;

  renderMarqueeBox = false; 

  SetZBufferParams();

  InitOpenGL();

  int windowWidth, windowHeight;
  GetClientSize(&windowWidth, &windowHeight);
  wxSizeEvent event(wxSize(windowWidth,windowHeight));
  OnSize(event);

  Refresh(TRUE);
}

MyGLCanvas::~MyGLCanvas()
{
  if (renderingMeshDisplayListAvailable)
  {
    // remove display list
    glDeleteLists(renderingMeshDisplayList, 1);
  }
  delete(m_context);
}

void MyGLCanvas::SetZBufferParams()
{
  zNear = cameraRadius * 0.02;
  zFar = cameraRadius * 20;
}

void MyGLCanvas::Reshape()
{
  int windowWidth, windowHeight;
  GetClientSize(&windowWidth, &windowHeight);
  SetZBufferParams();

/*
#ifndef __WXMOTIF__
  if ( GetContext() )
#endif
*/
  {
    glViewport(0, 0, windowWidth, windowHeight);

    glMatrixMode(GL_PROJECTION); // Select The Projection Matrix
    glLoadIdentity(); // Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, zNear, zFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }
}

void MyGLCanvas::OnPaint(wxPaintEvent& event)
{
  // must always be here
  //SetCurrent();
  wxGLCanvas::SetCurrent(*m_context);
  //wxPaintDC dc(this);
  wxPaintDC(this);

/*
  #ifndef __WXMOTIF__
    if (!GetContext()) 
      return;
  #endif
*/

  glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  glEnable(GL_DEPTH_TEST);

  if (uiState->vertexSelectionActivated)
  {
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  camera->Look(); // calls gluLookAt

  glDisable(GL_LIGHTING);

  if (renderingMeshDisplayListAvailable && showRenderingMesh)
    glCallList(renderingMeshDisplayList);

  if (precomputationState->linearModesAvailable && showLinearModes)
  {
    elapsedTime.StopCounter();
    time = elapsedTime.GetElapsedTime();
    double * ULinear = precomputationState->linearModalMatrix->GetMatrix();
    int n = precomputationState->linearModalMatrix->Getn();
    double omega = 0.5 * 2 * PI;
    for(int i=0; i<3*n; i++)
      uLinear[i] = linearRenderingMagnitude * sin(omega * time) * ULinear[ELT(3*n, i, renderedLinearMode)];

    RenderVolumetricMesh renderVolumetricMesh;

    if (uiState->showMaterialGroups)
      renderVolumetricMesh.SetDiscreteRenderingMode();
    else
      renderVolumetricMesh.SetFlatRenderingMode();

    renderVolumetricMesh.RenderSolidAndWireframeDeformation(precomputationState->simulationMesh, uLinear);
  }

  if (precomputationState->nonLinearModesAvailable && showNonLinearModes)
  {
    elapsedTime.StopCounter();
    time = elapsedTime.GetElapsedTime();
    double * UNonLinear = precomputationState->nonLinearModalMatrix->GetMatrix();
    int n = precomputationState->nonLinearModalMatrix->Getn();
    double omega = 0.5 * 2 * PI;
    for(int i=0; i<3*n; i++)
      uNonLinear[i] = nonLinearRenderingMagnitude * sin(omega * time) * UNonLinear[ELT(3*n, i, renderedNonLinearMode)];

    RenderVolumetricMesh renderVolumetricMesh;

    if (uiState->showMaterialGroups)
      renderVolumetricMesh.SetDiscreteRenderingMode();
    else
      renderVolumetricMesh.SetFlatRenderingMode();

    renderVolumetricMesh.RenderSolidAndWireframeDeformation(
      precomputationState->simulationMesh, uNonLinear);
  }

  if (uiState->vertexSelectionActivated)
  {
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    glStencilFunc(GL_ALWAYS, 1, ~0);
  }

  glColor3f(1.0,1.0,1.0);
  if (simulationMeshDisplayListAvailable && showSimulationMesh && !showLinearModes && uiState->renderMesh)
    glCallList(simulationMeshDisplayList);

  if (uiState->vertexSelectionActivated)
  {
    // render 0s into the stencil buffer from now on
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    glStencilFunc(GL_ALWAYS, 0, ~0);
  }

  // draw any extra 3D geometry here...
  // (none for now)

  // 2D primitives from this point onward

  if (uiState->vertexSelectionActivated)
  {
    // disable any stencil updates
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
    glStencilFunc(GL_ALWAYS, 1, ~0);
  }

  RenderVolumetricMesh renderVolumetricMesh;
  showFixedVertices = showSimulationMesh || showLinearModes || showNonLinearModes; // always show if showing the right mesh
  bool fixedVerticesAvailable = precomputationState->fixedVertices.size() > 0;
  if (fixedVerticesAvailable && showFixedVertices && (showSimulationMesh || showLinearModes || showNonLinearModes))
  {
    glPointSize(5.0);
    glColor3f(1,0,0);
    renderVolumetricMesh.RenderVertices(precomputationState->simulationMesh, 
      &(precomputationState->fixedVertices), false);
  }

  glLineWidth(2.0);

  if (uiState->showAxes)
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

  // on-screen rendering
  glDisable(GL_DEPTH_TEST);
  GLint view[4];
  glGetIntegerv (GL_VIEWPORT, view);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity(); 
  gluOrtho2D(0, view[2], 0, view[3]);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  if (renderMarqueeBox)
  {
    glLineWidth(2.0);
    glColor3f(0.5,0,0);
    glBegin(GL_LINE_LOOP);
    glVertex2f(marqueeBoxMin[0], marqueeBoxMin[1]);
    glVertex2f(marqueeBoxMax[0], marqueeBoxMin[1]);
    glVertex2f(marqueeBoxMax[0], marqueeBoxMax[1]);
    glVertex2f(marqueeBoxMin[0], marqueeBoxMax[1]);
    glEnd();
  }
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glFlush();
  SwapBuffers();
}

void MyGLCanvas::OnSize(wxSizeEvent& event)
{
  // this is necessary to update the context on some platforms
  //wxGLCanvas::OnSize(event);
  
  // Reset the OpenGL view aspect ratio
  Reshape();
}


void MyGLCanvas::OnEraseBackground(wxEraseEvent& event)
{
  // do nothing, to avoid flashing on MSW
}

void MyGLCanvas::OnMouse( wxMouseEvent& event )
{
  //wxSize sz(GetClientSize());
  if( event.ButtonDown() )
  {
    if (event.LeftDown())
    {
      GLint view[4];
      glGetIntegerv (GL_VIEWPORT, view);
      startDragX = event.GetX();
      startDragY = view[3]-1-event.GetY();
    }
    else if( event.MiddleIsDown() )
    {
    }
    else if( event.RightIsDown() )
    {
    }
  }
  else if (event.ButtonUp() )
  {
    if( event.LeftUp())
    {
      if (uiState->vertexSelectionActivated)
      {
        renderMarqueeBox = false;

        GLint view[4];
        glGetIntegerv (GL_VIEWPORT, view);

        int endDragX = event.GetX();
        int endDragY = view[3]-1-event.GetY();

        GLdouble model[16];
        glGetDoublev (GL_MODELVIEW_MATRIX, model);

        GLdouble proj[16];
        glGetDoublev (GL_PROJECTION_MATRIX, proj);

        if ( (endDragX - startDragX) * (endDragX - startDragX) + 
             (endDragY - startDragY) * (endDragY - startDragY) <= 2 )
        {
          // this was a single click
          float zValue;
          glReadPixels(startDragX,startDragY,1,1,GL_DEPTH_COMPONENT, GL_FLOAT,&zValue);

          GLubyte stencilValue;
          glReadPixels (startDragX, startDragY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, &stencilValue);

          GLdouble worldX, worldY, worldZ;
          gluUnProject (startDragX, startDragY, zValue, model, proj, view, &worldX, &worldY, &worldZ);

          if (stencilValue == 1)
          {
            int closestVertex = precomputationState->simulationMesh->getClosestVertex(Vec3d(worldX,worldY,worldZ));

            printf("Clicked on vertex: %d\n", closestVertex);

            if (event.ShiftDown())
            {
              if (precomputationState->fixedVertices.find(closestVertex) == precomputationState->fixedVertices.end())
                precomputationState->fixedVertices.insert(closestVertex);
              else
                precomputationState->fixedVertices.erase(closestVertex);
            }
            else if (event.ControlDown())
            {
              precomputationState->fixedVertices.insert(closestVertex);
            }
            else
            {
              precomputationState->fixedVertices.clear();
              precomputationState->fixedVertices.insert(closestVertex);
            }

            precomputationState->fixedVerticesAvailable = (precomputationState->fixedVertices.size() > 0);

            parent->UpdateMenus();

            Refresh(FALSE);
          }
          else
          {
          }
        }
        else
        {
          // this was a marquee-select drag

          // establish a sorted bounding box
          int bboxMin[2];
          int bboxMax[2];
          
          if (endDragX > startDragX)
          {
            bboxMin[0] = startDragX;
            bboxMax[0] = endDragX;
          }
          else
          {
            bboxMin[0] = endDragX;
            bboxMax[0] = startDragX;
          }

          if (endDragY > startDragY)
          {
            bboxMin[1] = startDragY;
            bboxMax[1] = endDragY;
          }
          else
          {
            bboxMin[1] = endDragY;
            bboxMax[1] = startDragY;
          }

          double cameraPosX, cameraPosY, cameraPosZ;
          camera->GetAbsWorldPosition(cameraPosX, cameraPosY, cameraPosZ);
          Vec3d cameraPos(cameraPosX, cameraPosY, cameraPosZ);

          GLdouble worldPos[3];
          Vec3d pyramid[4];
          int coords[4][2] = { 
            {bboxMin[0], bboxMin[1]},
            {bboxMax[0], bboxMin[1]},
            {bboxMax[0], bboxMax[1]},
            {bboxMin[0], bboxMax[1]}};

          for(int j=0; j<4; j++)
          {
            gluUnProject (coords[j][0], coords[j][1], 0.5, model, proj, view, 
              &worldPos[0], &worldPos[1], &worldPos[2]);
            pyramid[j][0] = worldPos[0] - cameraPos[0];
            pyramid[j][1] = worldPos[1] - cameraPos[1];
            pyramid[j][2] = worldPos[2] - cameraPos[2];
          }

          Vec3d faceNormal[4];
          for(int j=0; j<4; j++)
            faceNormal[j] = norm(cross(pyramid[(j+1)%4], pyramid[j]));

          if (!event.ShiftDown() && (!event.ControlDown()))
            precomputationState->fixedVertices.clear();

          for(int j=0; j<precomputationState->simulationMesh->getNumVertices(); j++)
          {
            Vec3d * vertexPos = precomputationState->simulationMesh->getVertex(j);
            Vec3d direction = (*vertexPos) - cameraPos;

            bool insidePyramid = true;
            for(int jj=0; jj<4; jj++)
              insidePyramid = insidePyramid && (dot(direction, faceNormal[jj]) >= 0);

            if (insidePyramid)
            {
              if (event.ShiftDown())
              {
                if (precomputationState->fixedVertices.find(j) == precomputationState->fixedVertices.end())
                  precomputationState->fixedVertices.insert(j);
                else
                  precomputationState->fixedVertices.erase(j);
              }
              else
              {
                precomputationState->fixedVertices.insert(j);
              }
            }
          } // end for all over vertices

          precomputationState->fixedVerticesAvailable = (precomputationState->fixedVertices.size() > 0);

          parent->UpdateMenus();

          Refresh(FALSE);
        } // end marquee selection drag
      }
    }
  }
  else if( event.Dragging() )
  {
    int mouseDeltaX = event.GetX() - mousePos[0];
    int mouseDeltaY = event.GetY() - mousePos[1];

    if (event.RightIsDown() )
    {
      camera->MoveRight(mouseDeltaX);
      camera->MoveUp(mouseDeltaY);
    }

    if (event.MiddleIsDown() || (event.LeftIsDown() && event.AltDown()) )
    {
      camera->ZoomIn(mouseDeltaY * cameraRadius * 0.05);
    }

    if (event.LeftIsDown() && (!event.AltDown()))
    {
      if (uiState->vertexSelectionActivated)
      {
        GLint view[4];
        glGetIntegerv (GL_VIEWPORT, view);

        int endDragX = event.GetX();
        int endDragY = view[3]-1-event.GetY();

        if ( (endDragX - startDragX) * (endDragX - startDragX) + 
             (endDragY - startDragY) * (endDragY - startDragY) > 2 )
        {
          renderMarqueeBox = true;  

          if (endDragX > startDragX)
          {
            marqueeBoxMin[0] = startDragX;
            marqueeBoxMax[0] = endDragX;
          }
          else
          {
            marqueeBoxMin[0] = endDragX;
            marqueeBoxMax[0] = startDragX;
          }

          if (endDragY > startDragY)
          {
            marqueeBoxMin[1] = startDragY;
            marqueeBoxMax[1] = endDragY;
          }
          else
          {
            marqueeBoxMin[1] = endDragY;
            marqueeBoxMax[1] = startDragY;
          }
        }
      }
    }

    Refresh(FALSE);
  }

  mousePos[0] = event.GetX();
  mousePos[1] = event.GetY();

  SetFocus();
}

void MyGLCanvas::OnKey(wxKeyEvent& event)
{
  switch( event.GetKeyCode() )
  {
    case WXK_RIGHT:
      camera->MoveFocusRight(-0.1 * camera->GetRadius());
      Refresh(FALSE);
    break;

    case WXK_LEFT:
      camera->MoveFocusRight(+0.1 * camera->GetRadius());
      Refresh(FALSE);
    break;

    case WXK_UP:
      camera->MoveFocusUp(-0.1 * camera->GetRadius());
      Refresh(FALSE);
    break;

    case WXK_DOWN:
      camera->MoveFocusUp(+0.1 * camera->GetRadius());
      Refresh(FALSE);
    break;

    case WXK_NUMPAD_ADD:
    break;

    case WXK_NUMPAD_SUBTRACT:
    break;

    case 'E':
      if (!event.ShiftDown())
      {
        uiState->renderMesh = !uiState->renderMesh;
        Refresh(FALSE);
      }
    break;

    case ']':
    break;

    default:
    break;
  }
}

void MyGLCanvas::InitOpenGL()
{
  // clear to white
  glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to light blue
  //glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  // clear to gray
  //glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);

  glEnable(GL_DEPTH_TEST);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  //glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
  //glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);

  //  glEnableClientState(GL_VERTEX_ARRAY);
  //  glEnableClientState(GL_NORMAL_ARRAY);

  //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Really Nice Perspective Calculations

  Reshape(); // calls gluPerspective, and other things
}

void MyGLCanvas::OnIdle(wxIdleEvent& event)
{
  if ((precomputationState->linearModesAvailable && showLinearModes) || 
      (precomputationState->nonLinearModesAvailable && showNonLinearModes))
  {
    this->Refresh(FALSE);
    event.RequestMore();
  }
}

void MyGLCanvas::DeleteRenderingMeshRenderData()
{
  if (renderingMeshDisplayListAvailable)
  {
    // remove old display list
    glDeleteLists(renderingMeshDisplayList, 1);
  }

  renderingMeshDisplayListAvailable = false;
}

void MyGLCanvas::UpdateRenderingMeshRenderData(int updateCamera)
{
  if (renderingMeshDisplayListAvailable)
  {
    // remove old display list
    glDeleteLists(renderingMeshDisplayList, 1);
  }

  renderingMeshDisplayListAvailable = true;

  // create display list
  renderingMeshDisplayList = glGenLists(1);
  glNewList(renderingMeshDisplayList, GL_COMPILE);
  {
    ObjMeshRender objMeshRender(precomputationState->renderingMesh);
    objMeshRender.setCustomColors(Vec3d(0,0,0.5));
    objMeshRender.render(OBJMESHRENDER_TRIANGLES, OBJMESHRENDER_CUSTOMCOLOR);
    glColor3f(0,0,0);

    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(3.0, 1000.0);
    objMeshRender.render(OBJMESHRENDER_EDGES, OBJMESHRENDER_NONE);
    glDisable(GL_POLYGON_OFFSET_LINE);
  }
  glEndList();

  if (updateCamera)
  {
    // set new camera parameters
    Vec3d bmin, bmax;
    precomputationState->renderingMesh->getBoundingBox(2.0, &bmin, &bmax);

    double center[3];
    center[0] = 0.5 * (bmin[0] + bmax[0]);
    center[1] = 0.5 * (bmin[1] + bmax[1]);
    center[2] = 0.5 * (bmin[2] + bmax[2]);

    cameraRadius = 0.6 * len(bmax-bmin);

    camera->SetFocusPosition(center);
    camera->SetPosition(cameraRadius, camera->GetPhi(), camera->GetTheta());

    Reshape();

    Refresh(false);
  }
}

void MyGLCanvas::UpdateSimulationMeshRenderData()
{
  if (simulationMeshDisplayListAvailable)
  {
    // remove old display list
    glDeleteLists(simulationMeshDisplayList, 1);
  }

  simulationMeshDisplayListAvailable = true;

  // create display list
  RenderVolumetricMesh renderVolumetricMesh;

  if (uiState->showMaterialGroups)
    renderVolumetricMesh.SetDiscreteRenderingMode();
  else
    renderVolumetricMesh.SetFlatRenderingMode();

  simulationMeshDisplayList = glGenLists(1);
  glNewList(simulationMeshDisplayList, GL_COMPILE);
    renderVolumetricMesh.RenderSolidAndWireframe(precomputationState->simulationMesh);
  glEndList();

  Vec3d center;
  double radius;
  precomputationState->simulationMesh->getMeshGeometricParameters(center, &radius);

  cameraRadius = 1.2 * 2 * radius;

  double centerv[3];
  center.convertToArray(centerv);
  camera->SetFocusPosition(centerv);

  camera->SetPosition(cameraRadius, camera->GetPhi(), camera->GetTheta());

  Reshape();

  Refresh(false);
}

void MyGLCanvas::UpdateLinearModesRenderData()
{
  time = 0.0;
  free(uLinear);
  uLinear = (double*) calloc (3 * precomputationState->linearModalMatrix->Getn(), sizeof(double));
}

void MyGLCanvas::UpdateNonLinearModesRenderData()
{
  time = 0.0;
  free(uNonLinear);
  uNonLinear = (double*) calloc (3 * precomputationState->nonLinearModalMatrix->Getn(), sizeof(double));
}

