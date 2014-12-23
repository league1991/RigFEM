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

// the functionality of the "View" menu

#include <wx/filename.h>

#include "largeModalDeformationFactory.h"

void MyFrame::UpdateModalToolbar()
{
  int newMode = modeSelectionControl->GetValue() - 1;

  if ((uiState.viewMode == UIState::VIEW_LINEAR_MODES) && precomputationState.linearModesAvailable)
  {
    myGLCanvas->SelectLinearMode(newMode);
  }

  if ((uiState.viewMode == UIState::VIEW_NONLINEAR_MODES) && precomputationState.nonLinearModesAvailable)
  {
    myGLCanvas->SelectNonLinearMode(newMode);
  }

  if (precomputationState.frequenciesAvailable && (uiState.viewMode == UIState::VIEW_LINEAR_MODES))
  {
    char s[4096];
    sprintf(s, "%f", precomputationState.frequencies[newMode]);
    frequencyTextCtrl->SetValue(wxString(s, wxConvUTF8));
  }
  else
  {
    frequencyTextCtrl->SetValue(_T("--"));
  }

  if ((uiState.viewMode == UIState::VIEW_LINEAR_MODES) && precomputationState.linearModesAvailable)
  {
    double renderingMagnitude;
    if ((modeAmplitudeControl->GetValue().ToDouble(&renderingMagnitude)))
      myGLCanvas->SetLinearRenderingMagnitude(renderingMagnitude);
  }
  else
    if ((uiState.viewMode == UIState::VIEW_NONLINEAR_MODES) && precomputationState.nonLinearModesAvailable)
    {
      double renderingMagnitude;
      if ((modeAmplitudeControl->GetValue().ToDouble(&renderingMagnitude)))
        myGLCanvas->SetNonLinearRenderingMagnitude(renderingMagnitude);
    }
    else
    {
      if (modeAmplitudeControl != NULL)
        modeAmplitudeControl->SetValue(_T("--"));
    }
}

void MyFrame::EnableModalToolbar(bool enable)
{
  modeSelectionControl->Enable(enable);
  modeSelectionText->Enable(enable);
  frequencyText->Enable(enable);
  frequencyTextCtrl->Enable(enable);
  modeAmplitudeText->Enable(enable);
  modeAmplitudeControl->Enable(enable);

  UpdateModalToolbar();
}

void MyFrame::ActivateVertexSelection(bool activate, int noViewSelect)
{
   uiState.vertexSelectionActivated = activate;
   menuBar->Check(ID_SelectFixedVertices, uiState.vertexSelectionActivated);

   if ((!noViewSelect) && (uiState.vertexSelectionActivated))
     SelectView(UIState::VIEW_SIMULATION_MESH);
}

void MyFrame::SelectView(UIState::viewModeType viewMode)
{
  myGLCanvas->ShowRenderingMesh(false);
  myGLCanvas->ShowSimulationMesh(false);
  myGLCanvas->ShowLinearModes(false);
  myGLCanvas->ShowNonLinearModes(false);

  this->uiState.viewMode = viewMode;
  switch(viewMode)
  {
    case UIState::VIEW_NONE:
      viewText->SetLabel(" View: No model loaded");
    break;

    case UIState::VIEW_RENDERING_MESH:
      myGLCanvas->ShowRenderingMesh(true);
      menuBar->Check(ID_ViewRenderingMesh, true);
      viewText->SetLabel(" View: Triangle mesh");
      EnableModalToolbar(false);
      ActivateVertexSelection(false);
    break;

    case UIState::VIEW_SIMULATION_MESH:
      myGLCanvas->ShowSimulationMesh(true);
      menuBar->Check(ID_ViewSimulationMesh, true);
      viewText->SetLabel(" View: Simulation mesh");
      ActivateVertexSelection(true, 1);
      EnableModalToolbar(false);
    break;

    case UIState::VIEW_LINEAR_MODES:
      myGLCanvas->ShowLinearModes(true);
      menuBar->Check(ID_ViewLinearModes, true);
      viewText->SetLabel(" View: Linear modes");
      EnableModalToolbar(true);
      modeSelectionControl->SetRange(1, precomputationState.rLin);
      if (modeSelectionControl->GetValue() > precomputationState.rLin)
        modeSelectionControl->SetValue(precomputationState.rLin);
      ActivateVertexSelection(false);
    break;

    case UIState::VIEW_NONLINEAR_MODES:
      myGLCanvas->ShowNonLinearModes(true);
      menuBar->Check(ID_ViewNonLinearModes, true);
      viewText->SetLabel(" View: Nonlinear modes");
      EnableModalToolbar(true);
      modeSelectionControl->SetRange(1, precomputationState.rNonLin);
      ActivateVertexSelection(false);
    break;

    default:
      ActivateVertexSelection(false);
    break;
  }

  myGLCanvas->Refresh();
}

void MyFrame::OnViewRenderingMesh(wxCommandEvent& event)
{
  SelectView(UIState::VIEW_RENDERING_MESH);
}

void MyFrame::OnViewSimulationMesh(wxCommandEvent& event)
{
  SelectView(UIState::VIEW_SIMULATION_MESH);
}

void MyFrame::OnViewLinearModes(wxCommandEvent& event)
{
  SelectView(UIState::VIEW_LINEAR_MODES);
}

void MyFrame::OnViewNonlinearModes(wxCommandEvent& event)
{
  SelectView(UIState::VIEW_NONLINEAR_MODES);
}

void MyFrame::OnViewRuntimeSimulation(wxCommandEvent& event)
{
  SelectView(UIState::VIEW_RUNTIME_SIMULATION);
}

void MyFrame::OnShowAxes(wxCommandEvent& event)
{
  uiState.showAxes = event.IsChecked();
  myGLCanvas->Refresh(FALSE);
}

void MyFrame::OnShowMaterialGroups(wxCommandEvent& event)
{
  uiState.showMaterialGroups = event.IsChecked();
  if (precomputationState.simulationMeshAvailable)
    myGLCanvas->UpdateSimulationMeshRenderData();
  myGLCanvas->Refresh(FALSE);
}

void MyFrame::UpdateMenus()
{
  menuBar->Enable(ID_ViewRenderingMesh, precomputationState.renderingMeshAvailable);
  menuBar->Enable(ID_ViewSimulationMesh, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ViewLinearModes, precomputationState.linearModesAvailable);
  menuBar->Enable(ID_ViewNonLinearModes, precomputationState.nonLinearModesAvailable);

  menuBar->Enable(ID_LoadMesh, true);
  menuBar->Enable(ID_RemoveTriangleMesh, precomputationState.renderingMeshAvailable);
  menuBar->Enable(ID_LoadVoxelization, true);
  menuBar->Enable(ID_InterpolateLinearModes, 
    precomputationState.renderingMeshAvailable && precomputationState.simulationMeshAvailable && precomputationState.linearModesAvailable);
  menuBar->Enable(ID_InterpolateNonLinearModes, 
    precomputationState.renderingMeshAvailable && precomputationState.simulationMeshAvailable && precomputationState.nonLinearModesAvailable);
  menuBar->Enable(ID_InterpolateModalDerivatives, 
    precomputationState.renderingMeshAvailable && precomputationState.simulationMeshAvailable && precomputationState.modalDerivativesAvailable);
  menuBar->Enable(ID_InterpolateSimulationData, 
    precomputationState.renderingMeshAvailable && precomputationState.simulationMeshAvailable && precomputationState.sketchDataAvailable);
  menuBar->Enable(ID_SaveInterpolant, 
    precomputationState.renderingMeshAvailable && precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_Voxelize, precomputationState.renderingMeshAvailable);

  menuBar->Enable(ID_LoadFixedVertices, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_SaveVoxelization, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ExportSurfaceMesh, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_MaterialProperties, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_SelectFixedVertices, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ClearFixedVertices, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_FixedVerticesInformation, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_LoadLinearModes, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_LoadFrequencies, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_MeshInformation, precomputationState.renderingMeshAvailable || precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_SaveFixedVertices, precomputationState.fixedVerticesAvailable);
  menuBar->Enable(ID_ComputeLinearModes, precomputationState.simulationMeshAvailable );
  menuBar->Enable(ID_SaveLinearModes, precomputationState.linearModesAvailable);

  menuBar->Enable(ID_ExportMassMatrixFull, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ExportStiffnessMatrixFull, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ExportMassMatrixConstrained, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_ExportStiffnessMatrixConstrained, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_DisplayFrequencies, precomputationState.frequenciesAvailable);
  menuBar->Enable(ID_SaveFrequencies, precomputationState.frequenciesAvailable);
  menuBar->Enable(ID_ScaleTo1Hz, precomputationState.frequenciesAvailable);
  menuBar->Enable(ID_ScaleFrequencies, precomputationState.frequenciesAvailable);
  menuBar->Enable(ID_DeleteFrequencies, precomputationState.linearModesAvailable && precomputationState.frequenciesAvailable);

  menuBar->Enable(ID_DisplayRigidModes, precomputationState.linearModesAvailable && precomputationState.frequenciesAvailable);
  menuBar->Enable(ID_SetRigidModes, precomputationState.linearModesAvailable && precomputationState.frequenciesAvailable);

  menuBar->Enable(ID_LoadNonLinearModes, precomputationState.simulationMeshAvailable);
  menuBar->Enable(ID_SaveNonLinearModes, precomputationState.nonLinearModesAvailable);

  menuBar->Enable(ID_ComputeNonLinearModes, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_LoadModalDerivatives, precomputationState.linearModesAvailable);
  menuBar->Enable(ID_SaveModalDerivatives, precomputationState.modalDerivativesAvailable);

  menuBar->Enable(ID_ComputeModalDerivatives, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_LoadSketchData, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_LoadCubicPolynomials, precomputationState.nonLinearModesAvailable);
  menuBar->Enable(ID_SaveCubicPolynomials, precomputationState.cubicPolynomialsAvailable);

  menuBar->Enable(ID_ComputeCubicPolynomials, precomputationState.simulationMeshAvailable);

  menuBar->Enable(ID_BatchPreprocessSingleSimulation, precomputationState.simulationMeshAvailable );
  menuBar->Enable(ID_BatchPreprocessManySimulations, true);

  menuBar->Enable(ID_PrepareRuntimeSimulation, precomputationState.simulationMeshAvailable && precomputationState.nonLinearModesAvailable && precomputationState.cubicPolynomialsAvailable);
  menuBar->Enable(ID_TuneRuntimeSimulation, precomputationState.simulationMeshAvailable && precomputationState.nonLinearModesAvailable && precomputationState.cubicPolynomialsAvailable && precomputationState.runtimeSimReady);
  menuBar->Enable(ID_LaunchRuntimeSimulation, precomputationState.simulationMeshAvailable && precomputationState.nonLinearModesAvailable && precomputationState.cubicPolynomialsAvailable && precomputationState.runtimeSimReady);
}

void MyFrame::SaveCurrentWorkingDirectory(wxString & completePathIncludingFilename)
{
  wxString volume, path, name, ext;
  wxFileName::SplitPath(completePathIncludingFilename, &volume, &path, &name, &ext);

  #ifdef WIN32
    uiState.currentWorkingDirectory = volume + _T(":") + path;
  #else
    uiState.currentWorkingDirectory = volume + path;
  #endif

  FILE * fout = fopen("directory.log", "w");
  fprintf(fout, "%s\n", (const char*)uiState.currentWorkingDirectory.mb_str());
  fclose(fout);
}

