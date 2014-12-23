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

#ifndef _CANVAS_H_
#define _CANVAS_H_

#include "wx/glcanvas.h"
#include "camera.h"
#include "objMesh.h"
#include "cubicMesh.h"
#include "modalMatrix.h"
#include "performanceCounter.h"
#include "largeModalDeformationFactory.h"
#include "states.h"

// this class does all the OpenGL rendering in this program
class MyGLCanvas : public wxGLCanvas
{
public:

  MyGLCanvas(PrecomputationState * precomputationState, UIState * uiState,
    MyFrame * parent, wxWindowID id = -1, 
    const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, 
    long style=0, const wxString& name="GLCanvas", int* attribList = 0, 
    const wxPalette& palette = wxNullPalette);
  ~MyGLCanvas();

  void OnPaint(wxPaintEvent& event);
  void OnSize(wxSizeEvent& event);
  void OnIdle(wxIdleEvent& event);
  void OnKey(wxKeyEvent& event);
  void OnMouse( wxMouseEvent& event);
  void OnEraseBackground(wxEraseEvent& event);

  void InitOpenGL();
  void Reshape();
  void UpdateRenderingMeshRenderData(int updateCamera=1);
  void DeleteRenderingMeshRenderData();
  void UpdateSimulationMeshRenderData();
  void UpdateLinearModesRenderData();
  void UpdateNonLinearModesRenderData();

  inline SphericalCamera * GetCamera() { return camera; }

  void SelectLinearMode(int mode) { renderedLinearMode = mode; }
  void SelectNonLinearMode(int mode) { renderedNonLinearMode = mode; }
  void SetLinearRenderingMagnitude(double linearRenderingMagnitude) { this->linearRenderingMagnitude = linearRenderingMagnitude; }
  void SetNonLinearRenderingMagnitude(double nonLinearRenderingMagnitude) 
    { this->nonLinearRenderingMagnitude = nonLinearRenderingMagnitude; }

  // switches to enable/disable display of certain components
  inline void ShowRenderingMesh( bool showRenderingMesh ) { this->showRenderingMesh = showRenderingMesh;}
  inline void ShowSimulationMesh( bool showSimulationMesh ) { this->showSimulationMesh = showSimulationMesh;}
  inline void ShowLinearModes( bool showLinearModes ) { this->showLinearModes = showLinearModes;}
  inline void ShowNonLinearModes( bool showNonLinearModes ) { this->showNonLinearModes = showNonLinearModes;}

protected:
  MyFrame * parent;
  PrecomputationState * precomputationState;
  UIState * uiState;

  wxGLContext*	m_context;

  double zNear, zFar;
  double cameraRadius;
  SphericalCamera * camera;
  void SetZBufferParams();

  double time;
  PerformanceCounter elapsedTime;

  double mousePos[2];
  int startDragX, startDragY;
  bool renderMarqueeBox;  
  int marqueeBoxMin[2];
  int marqueeBoxMax[2];

  bool renderingMeshDisplayListAvailable;
  bool showRenderingMesh;
  GLuint renderingMeshDisplayList;

  bool simulationMeshDisplayListAvailable;
  GLuint simulationMeshDisplayList;
  bool showSimulationMesh;

  bool showFixedVertices;

  bool showLinearModes;
  double linearRenderingMagnitude;
  double * uLinear;

  bool showNonLinearModes;
  double nonLinearRenderingMagnitude;
  double * uNonLinear;

  int renderedLinearMode;
  int renderedNonLinearMode;

  DECLARE_EVENT_TABLE()
};

#endif

