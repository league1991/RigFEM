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

/*
  See instructions.txt for instructions on how to use this program.
*/

#ifndef _LARGEMODALDEFORMATIONFACTORY_H_
#define _LARGEMODALDEFORMATIONFACTORY_H_

#include <set>
using namespace std;

#include "wx/wx.h" 
#include "wx/spinctrl.h" 

#include "objMesh.h"
#include "loadList.h"
#include "volumetricMesh.h"
#include "modalMatrix.h"
#include "StVKReducedInternalForces.h"

#include "states.h"
class MyGLCanvas;
class MyFrame;
#include "canvas.h"

class MyApp: public wxApp
{
  virtual bool OnInit();
};

// the main class in this program
class MyFrame: public wxFrame
{
public:

  MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

  void OnMenu(wxCommandEvent& event);
  void OnMenuClose(wxMenuEvent& event);
  void OnMenuOpen(wxMenuEvent& event);

  void OnLoadProject(wxCommandEvent& event);
  void OnSaveProject(wxCommandEvent& event);
  void OnQuit(wxCommandEvent& event);

  void OnViewRenderingMesh(wxCommandEvent& event);
  void OnViewSimulationMesh(wxCommandEvent& event);
  void OnViewLinearModes(wxCommandEvent& event);
  void OnViewNonlinearModes(wxCommandEvent& event);
  void OnViewRuntimeSimulation(wxCommandEvent& event);
  void OnShowAxes(wxCommandEvent& event);
  void OnShowMaterialGroups(wxCommandEvent& event);

  void OnLoadMesh(wxCommandEvent& event);
  void OnRemoveTriangleMesh(wxCommandEvent& event);
  void OnVoxelize(wxCommandEvent& event);

  void OnInterpolateLinearModes(wxCommandEvent& event);
  void OnInterpolateNonLinearModes(wxCommandEvent& event);
  void OnInterpolateModalDerivatives(wxCommandEvent& event);
  void OnInterpolateSimulationData(wxCommandEvent& event);
  void OnSaveInterpolant(wxCommandEvent& event);

  void OnLoadSimulationMesh(wxCommandEvent& event);
  void OnSaveSimulationMesh(wxCommandEvent& event);
  void OnExportSurfaceMesh(wxCommandEvent& event);
  void OnMaterialProperties(wxCommandEvent& event);
  void OnMeshInformation(wxCommandEvent& event);

  void OnLoadFixedVertices(wxCommandEvent& event);
  void OnSaveFixedVertices(wxCommandEvent& event);
  void OnSelectFixedVertices(wxCommandEvent& event);
  void OnClearFixedVertices(wxCommandEvent& event);
  void OnFixedVerticesInformation(wxCommandEvent& event);
 
  void OnLoadLinearModes(wxCommandEvent& event);
  void OnSaveLinearModes(wxCommandEvent& event);
  void OnExportLinearModes(wxCommandEvent& event);
  void OnComputeLinearModes(wxCommandEvent& event);
  void OnExportMassMatrixFull(wxCommandEvent& event);
  void OnExportStiffnessMatrixFull(wxCommandEvent& event);
  void OnExportMassMatrixConstrained(wxCommandEvent& event);
  void OnExportStiffnessMatrixConstrained(wxCommandEvent& event);

  void OnLoadFrequencies(wxCommandEvent& event);
  void OnSaveFrequencies(wxCommandEvent& event);
  void OnExportFrequencies(wxCommandEvent& event);
  void OnDisplayFrequencies(wxCommandEvent& event);
  void OnScaleTo1Hz(wxCommandEvent& event);
  void OnScaleFrequencies(wxCommandEvent& event);
  void OnDeleteFrequencies(wxCommandEvent&);
  void OnDisplayRigidModes(wxCommandEvent&);
  void OnSetRigidModes(wxCommandEvent&);

  void OnLoadNonLinearModes(wxCommandEvent& event);
  void OnSaveNonLinearModes(wxCommandEvent& event);
  void OnExportNonLinearModes(wxCommandEvent& event);
  void OnComputeNonLinearModes(wxCommandEvent& event);

  void OnLoadSketchData(wxCommandEvent & event);

  void OnLoadModalDerivatives(wxCommandEvent& event);
  void OnSaveModalDerivatives(wxCommandEvent& event);
  void OnExportModalDerivatives(wxCommandEvent& event);
  void OnComputeModalDerivatives(wxCommandEvent& event);

  void OnLoadCubicPolynomials(wxCommandEvent& event);
  void OnSaveCubicPolynomials(wxCommandEvent& event);
  void OnComputeCubicPolynomials(wxCommandEvent& event);

  void OnBatchPreprocessSingleSimulation(wxCommandEvent& event);
  void OnBatchPreprocessManySimulations(wxCommandEvent& event);

  void OnPrepareRuntimeSimulation(wxCommandEvent& event);
  void OnLaunchRuntimeSimulation(wxCommandEvent& event);
  void OnTuneRuntimeSimulation(wxCommandEvent& event);

  void OnConvertTextMatrixToBinaryMatrix(wxCommandEvent& event);
  void OnConvertBinaryMatrixToTextMatrix(wxCommandEvent& event);

  void OnAbout(wxCommandEvent& event);
  void OnMouseHelp(wxCommandEvent& event);
  void OnExampleHelp(wxCommandEvent& event);

  void OnChangeRenderedMode(wxSpinEvent& event);
  void OnChangeRenderingAmplitude(wxCommandEvent& event);

  DECLARE_EVENT_TABLE()

protected:

  wxMenuBar * menuBar;
  wxToolBar * toolBar;
  MyGLCanvas * myGLCanvas;

  wxStaticText * viewText;
  wxSpinCtrl * modeSelectionControl;
  wxStaticText * modeSelectionText;
  wxStaticText * frequencyText;
  wxTextCtrl * frequencyTextCtrl;
  wxStaticText * modeAmplitudeText;
  wxTextCtrl * modeAmplitudeControl;
  wxTextCtrl * inputMatrixTextCtrl;
  wxTextCtrl * outputMatrixTextCtrl;
  wxString versionString;

  PrecomputationState precomputationState; // the current state of the precomputation
  UIState uiState; // the current state of the user interface

  void UpdateMenus();
  void CreateOpenGLWindow();
  void InitOpenGL();
  string int2string(int n);
  string double2string(double d);
  double MaxModeMagnitude(ModalMatrix * modalMatrix);
  void DeleteRenderingMesh();
  int LoadSimulationMesh(wxString & meshFilename);
  void ComputeSimulationMeshGeometricParameters();
  void SetAutoRenderingMagnitude(ModalMatrix * modalMatrix);
  void SetDefaultMaterialParameters();
  void CreateRenderingMeshFromSimulationMesh();
  void ActivateVertexSelection(bool activate, int noViewSelect=0);
  int LoadFixedVertices(wxString & fixedVerticesFilename);
  void RemoveSixRigidModes(int numVectors, double * x);
  void ScaleYoungsModulus(double factor);
  void ExportMassMatrix(bool fullMatrix);
  void ExportStiffnessMatrix(bool fullMatrix);
  void errMsg(wxString title, wxString message, wxWindow* parent = NULL );

  // viewing
  void SelectView(UIState::viewModeType viewMode);
  void EnableModalToolbar(bool enable);
  void UpdateModalToolbar();

  // interpolation
  void BuildInterpolant();
  void InterpolateMatrix(wxString dataDescription, ModalMatrix * inputModalMatrix);
  void DeleteInterpolationData();
  void DeallocateSimulationData();

  // "worker" routines
  void * VoxelizationWorker(ObjMesh * objMesh, int resolution, bool fillHoles);
  void * LinearModesWorker(
      int numDesiredModes,
      int * r, double ** frequencies, double ** linearModes );
  void * NonLinearModesWorker(
      int * code, int dataOrigin, int numNonLinearModes, double ** modes_ );
  void ComputeModalDerivatives(int * code, double ** modalDerivatives);
  void * CubicPolynomialsWorker(int * code, StVKReducedInternalForces ** newCubicPolynomial);

  // runtime
  int PreprocessSingleSimulation(wxString & projectName, int numLinearModes, int numNonLinearModes, int numComputationThreads);
  void PreprocessManySimulationsDeallocateHelper(int numModels, char ** meshFiles, char ** bouFiles, int * numLinearModes, int * numNonLinearModes);

  // threads (not used) (in the future, this could be used to add the ability to i nterrupt long computations)
  int computationRunning;
  int comparator(const void * a, const void * b);
  void SpawnComputationThread(void * (MyFrame::*computationThread) (void*));
  void StopComputation(int parameter);
  void * VoxelizationThread(void * vargp);

  void SaveCurrentWorkingDirectory(wxString & completePathIncludingFilename);

  friend class MyGLCanvas;
  friend class PreprocessManySimulationsDialog;
};

enum
{
    ID_Quit = 1,
    ID_LoadProject,
    ID_SaveProject,
    ID_About,
    ID_MouseHelp,
    ID_ExampleHelp,

    ID_ViewRenderingMesh,
    ID_ViewSimulationMesh,
    ID_ViewLinearModes,
    ID_ViewNonLinearModes,
    ID_ShowAxes,
    ID_ShowMaterialGroups,

    ID_LoadMesh,
    ID_RemoveTriangleMesh,
    ID_Voxelize,

    ID_InterpolateLinearModes, 
    ID_InterpolateNonLinearModes, 
    ID_InterpolateModalDerivatives,
    ID_InterpolateSimulationData,
    ID_SaveInterpolant,

    ID_LoadVoxelization,
    ID_SaveVoxelization,
    ID_ExportSurfaceMesh,
    ID_MaterialProperties,
    ID_MeshInformation,

    ID_LoadFixedVertices,
    ID_SaveFixedVertices,
    ID_SelectFixedVertices,
    ID_ClearFixedVertices,
    ID_FixedVerticesInformation,

    ID_LoadLinearModes,
    ID_SaveLinearModes,
    ID_ExportLinearModes,
    ID_ComputeLinearModes,
    ID_ExportMassMatrixFull,
    ID_ExportMassMatrixConstrained,
    ID_ExportStiffnessMatrixFull,
    ID_ExportStiffnessMatrixConstrained,
    ID_LoadFrequencies,
    ID_SaveFrequencies,
    ID_ExportFrequencies,
    ID_DisplayFrequencies,
    ID_ScaleTo1Hz,
    ID_ScaleFrequencies,
    ID_DeleteFrequencies,
    ID_DisplayRigidModes,
    ID_SetRigidModes,

    ID_LoadNonLinearModes,
    ID_SaveNonLinearModes,
    ID_ExportNonLinearModes,
    ID_ComputeNonLinearModes,
    ID_LoadModalDerivatives,
    ID_SaveModalDerivatives,
    ID_ExportModalDerivatives,
    ID_ComputeModalDerivatives,
    ID_LoadSketchData,

    ID_LoadCubicPolynomials,
    ID_SaveCubicPolynomials,
    ID_ComputeCubicPolynomials,
    ID_BatchPreprocessSingleSimulation,
    ID_BatchPreprocessManySimulations,
    ID_ScriptFilenameButton,
    ID_PrepareRuntimeSimulation,
    ID_LaunchRuntimeSimulation,
    ID_TuneRuntimeSimulation,

    ID_ConvertTextMatrixToBinaryMatrix,
    ID_ConvertBinaryMatrixToTextMatrix,

    ID_ModeSelection,
    ID_RenderingAmplitude,

    ID_inputMatrixFilenameButton,
    ID_outputMatrixFilenameButton,

    ID_MenuLast
};

// a dialog wrapper
template <typename Dialog>
class WxDialogWrapper
{
public:
  WxDialogWrapper( Dialog* dialog = 0 ): dlg(dialog) {}
  ~WxDialogWrapper() { reset(); }

  void reset(Dialog* dlg_ = 0) { if (dlg) dlg->Destroy(); dlg = dlg_; }

  Dialog* operator->() { return dlg; }
  Dialog& operator*() { return *dlg; }

  const Dialog* get() const { return dlg; }
  Dialog* get() { return dlg; }

protected:
  Dialog * dlg;
};

#endif

