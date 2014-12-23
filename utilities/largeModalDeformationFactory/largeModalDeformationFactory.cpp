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

// initialization routines

#include "largeModalDeformationFactory.h"
#include "icon.cpp"
#include <wx/filename.h>
using namespace std;

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
  EVT_MENU_OPEN(MyFrame::OnMenuOpen)
  EVT_MENU_CLOSE(MyFrame::OnMenuClose)

  EVT_MENU(ID_Quit, MyFrame::OnQuit)
  EVT_MENU(ID_LoadProject, MyFrame::OnLoadProject)
  EVT_MENU(ID_SaveProject, MyFrame::OnSaveProject)

  EVT_MENU(ID_ViewRenderingMesh, MyFrame::OnViewRenderingMesh)
  EVT_MENU(ID_ViewSimulationMesh, MyFrame::OnViewSimulationMesh)
  EVT_MENU(ID_ViewLinearModes, MyFrame::OnViewLinearModes)
  EVT_MENU(ID_ViewNonLinearModes, MyFrame::OnViewNonlinearModes)
  EVT_MENU(ID_ShowAxes, MyFrame::OnShowAxes)
  EVT_MENU(ID_ShowMaterialGroups, MyFrame::OnShowMaterialGroups)

  EVT_MENU(ID_LoadMesh, MyFrame::OnLoadMesh)
  EVT_MENU(ID_RemoveTriangleMesh, MyFrame::OnRemoveTriangleMesh)
  EVT_MENU(ID_Voxelize, MyFrame::OnVoxelize)

  EVT_MENU(ID_InterpolateLinearModes, MyFrame::OnInterpolateLinearModes)
  EVT_MENU(ID_InterpolateNonLinearModes, MyFrame::OnInterpolateNonLinearModes)
  EVT_MENU(ID_InterpolateModalDerivatives, MyFrame::OnInterpolateModalDerivatives)
  EVT_MENU(ID_InterpolateSimulationData, MyFrame::OnInterpolateSimulationData)
  EVT_MENU(ID_SaveInterpolant, MyFrame::OnSaveInterpolant)

  EVT_MENU(ID_LoadVoxelization, MyFrame::OnLoadSimulationMesh)
  EVT_MENU(ID_SaveVoxelization, MyFrame::OnSaveSimulationMesh)
  EVT_MENU(ID_ExportSurfaceMesh, MyFrame::OnExportSurfaceMesh)
  EVT_MENU(ID_MaterialProperties, MyFrame::OnMaterialProperties)
  EVT_MENU(ID_MeshInformation, MyFrame::OnMeshInformation)

  EVT_MENU(ID_LoadFixedVertices, MyFrame::OnLoadFixedVertices)
  EVT_MENU(ID_SaveFixedVertices, MyFrame::OnSaveFixedVertices)
  EVT_MENU(ID_SelectFixedVertices, MyFrame::OnSelectFixedVertices)
  EVT_MENU(ID_ClearFixedVertices, MyFrame::OnClearFixedVertices)
  EVT_MENU(ID_FixedVerticesInformation, MyFrame::OnFixedVerticesInformation)

  EVT_MENU(ID_LoadLinearModes, MyFrame::OnLoadLinearModes)
  EVT_MENU(ID_SaveLinearModes, MyFrame::OnSaveLinearModes)
  EVT_MENU(ID_ComputeLinearModes, MyFrame::OnComputeLinearModes)
  EVT_MENU(ID_ExportMassMatrixFull, MyFrame::OnExportMassMatrixFull)
  EVT_MENU(ID_ExportStiffnessMatrixFull, MyFrame::OnExportStiffnessMatrixFull)
  EVT_MENU(ID_ExportMassMatrixConstrained, MyFrame::OnExportMassMatrixConstrained)
  EVT_MENU(ID_ExportStiffnessMatrixConstrained, MyFrame::OnExportStiffnessMatrixConstrained)

  EVT_MENU(ID_LoadFrequencies, MyFrame::OnLoadFrequencies)
  EVT_MENU(ID_SaveFrequencies, MyFrame::OnSaveFrequencies)
  EVT_MENU(ID_ExportFrequencies, MyFrame::OnExportFrequencies)
  EVT_MENU(ID_DisplayFrequencies, MyFrame::OnDisplayFrequencies)
  EVT_MENU(ID_ScaleTo1Hz, MyFrame::OnScaleTo1Hz)
  EVT_MENU(ID_ScaleFrequencies, MyFrame::OnScaleFrequencies)
  EVT_MENU(ID_DeleteFrequencies, MyFrame::OnDeleteFrequencies)

  EVT_MENU(ID_DisplayRigidModes, MyFrame::OnDisplayRigidModes)
  EVT_MENU(ID_SetRigidModes, MyFrame::OnSetRigidModes)

  EVT_MENU(ID_LoadNonLinearModes, MyFrame::OnLoadNonLinearModes)
  EVT_MENU(ID_SaveNonLinearModes, MyFrame::OnSaveNonLinearModes)
  EVT_MENU(ID_ExportNonLinearModes, MyFrame::OnExportNonLinearModes)
  EVT_MENU(ID_ComputeNonLinearModes, MyFrame::OnComputeNonLinearModes)

  EVT_MENU(ID_LoadSketchData,  MyFrame::OnLoadSketchData)

  EVT_MENU(ID_LoadModalDerivatives,  MyFrame::OnLoadModalDerivatives)
  EVT_MENU(ID_SaveModalDerivatives,  MyFrame::OnSaveModalDerivatives)
  EVT_MENU(ID_ExportModalDerivatives,  MyFrame::OnExportModalDerivatives)
  EVT_MENU(ID_ComputeModalDerivatives,  MyFrame::OnComputeModalDerivatives)

  EVT_MENU(ID_LoadCubicPolynomials, MyFrame::OnLoadCubicPolynomials)
  EVT_MENU(ID_SaveCubicPolynomials, MyFrame::OnSaveCubicPolynomials)
  EVT_MENU(ID_ComputeCubicPolynomials, MyFrame::OnComputeCubicPolynomials)

  EVT_MENU(ID_BatchPreprocessSingleSimulation, MyFrame::OnBatchPreprocessSingleSimulation)
  EVT_MENU(ID_BatchPreprocessManySimulations, MyFrame::OnBatchPreprocessManySimulations)

  EVT_MENU(ID_PrepareRuntimeSimulation, MyFrame::OnPrepareRuntimeSimulation)
  EVT_MENU(ID_TuneRuntimeSimulation, MyFrame::OnTuneRuntimeSimulation)
  EVT_MENU(ID_LaunchRuntimeSimulation, MyFrame::OnLaunchRuntimeSimulation)

  EVT_MENU(ID_ConvertTextMatrixToBinaryMatrix, MyFrame::OnConvertTextMatrixToBinaryMatrix)
  EVT_MENU(ID_ConvertBinaryMatrixToTextMatrix, MyFrame::OnConvertBinaryMatrixToTextMatrix)

  EVT_MENU(ID_About, MyFrame::OnAbout)
  EVT_MENU(ID_MouseHelp, MyFrame::OnMouseHelp)
  EVT_MENU(ID_ExampleHelp, MyFrame::OnExampleHelp)

  EVT_SPINCTRL(ID_ModeSelection, MyFrame::OnChangeRenderedMode)
  EVT_TEXT(ID_RenderingAmplitude, MyFrame::OnChangeRenderingAmplitude)

END_EVENT_TABLE()

bool MyApp::OnInit()
{
  //#ifndef __WXMSW__
  #ifdef __APPLE__
    if (freopen("stdout.txt", "w", stdout) == NULL)
      printf("Warning: failed to redirect output to stdout.txt\n");
    if (freopen("stderr.txt", "w", stderr) == NULL)
      printf("Warning: failed to redirect output to stderr.txt\n");
  #endif

  MyFrame *frame = new MyFrame( _T("Large Modal Deformation Factory"), 
    wxPoint(50,50), wxSize(800,600) );
  frame->Show(TRUE);
  SetTopWindow(frame);
  return TRUE;
} 

MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
: wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
  versionString = _T("3.4");
  printf("Large modal deformation factory. Version %s .\n", (const char*)versionString.mb_str());
  printf("By Jernej Barbic, CMU, MIT, USC, 2007-2014.\n");

  printf("Initializing application...");
  precomputationState.renderingMeshAvailable = false;
  precomputationState.simulationMeshAvailable = false;
  precomputationState.fixedVerticesAvailable = false;
  precomputationState.linearModesAvailable = false;
  precomputationState.frequenciesAvailable = false;
  precomputationState.nonLinearModesAvailable = false;
  precomputationState.cubicPolynomialsAvailable = false;

  precomputationState.renderingMesh = NULL;

  uiState.cubicMeshResolution = 65;
  uiState.fillInteriorChambers = false;
  uiState.showMaterialGroups = true;
  precomputationState.simulationMesh = NULL;
  precomputationState.FFDinterpolationVertices = NULL;
  precomputationState.FFDinterpolationWeights = NULL;
  precomputationState.simulationSurfaceMesh = NULL;
  precomputationState.fixedVertices.empty();

  uiState.projectName = _T("defaultProject");
  uiState.renderingMeshFilename = _T("none");

  precomputationState.linearModalMatrix = NULL;
  precomputationState.rLin = 0;
  uiState.numComputedLinearModes = 10;
  uiState.firstModalComputation = true;

  precomputationState.frequencies = NULL;
  uiState.lastFrequencyScaleFactor = 1.0;
  uiState.eraseRangeLo = 0;
  uiState.eraseRangeHi = 10;
  precomputationState.numRigidModes = 0;

  precomputationState.nonLinearModalMatrix = NULL;
  precomputationState.rNonLin = 0;
  uiState.numComputedNonLinearModes = 20;

  precomputationState.modalDerivativesAvailable = false;
  precomputationState.numDeriv = 0;
  precomputationState.modalDerivativesMatrix = NULL;

  precomputationState.sketchDataAvailable = false;
  precomputationState.sketchDataMatrix = NULL;

  precomputationState.cubicPolynomials = NULL;
  uiState.numComputationThreads = wxThread::GetCPUCount();
  if (uiState.numComputationThreads <= 0)
    uiState.numComputationThreads = 1;

  precomputationState.runtimeSimReady = false;

  precomputationState.interpolationDataAvailable = false;
  precomputationState.interpolationData_vertices = NULL;
  precomputationState.interpolationData_weights = NULL;

  uiState.showAxes = true;
  uiState.viewMode = UIState::VIEW_NONE;

  uiState.manySimulationsscriptFilename = _T("none");

  uiState.scaleRealTimeTo1HzCheckBox = true;
  
  if (wxFileName::FileExists("directory.log"))
  {
    FILE * fin = fopen("directory.log", "r");
    char directoryName[4096];
    int code = fscanf(fin, "%s", directoryName);
    if (code != 0)
      uiState.currentWorkingDirectory = wxString(directoryName);
    else
      uiState.currentWorkingDirectory = _T("");
    fclose(fin);
  }
  else
    uiState.currentWorkingDirectory = _T("");

  modeAmplitudeControl = NULL;

  wxMenu * menuFile = new wxMenu;
  menuFile->Append( ID_Quit, _T("E&xit") );

  wxMenu * menuView = new wxMenu;
  menuView->AppendRadioItem(ID_ViewRenderingMesh, _T("&Triangle mesh"));
  menuView->AppendRadioItem(ID_ViewSimulationMesh, _T("&Simulation mesh"));
  menuView->AppendRadioItem(ID_ViewLinearModes, _T("&Linear modes"));
  menuView->AppendRadioItem(ID_ViewNonLinearModes, _T("&Nonlinear modes"));
  menuView->AppendSeparator();
  menuView->AppendCheckItem( ID_ShowAxes, _T("& Coordinate axes"));
  menuView->AppendCheckItem( ID_ShowMaterialGroups, _T("& Material groups"));

  wxMenu * menuFixedVertices = new wxMenu;
  menuFixedVertices->Append( ID_LoadFixedVertices, _T("&Load...") );
  menuFixedVertices->Append( ID_SaveFixedVertices, _T("&Save...") );
  menuFixedVertices->AppendSeparator();
  menuFixedVertices->AppendCheckItem( ID_SelectFixedVertices, _T("&Select") );
  menuFixedVertices->Append( ID_ClearFixedVertices, _T("&Clear") );
  menuFixedVertices->AppendSeparator();
  menuFixedVertices->Append( ID_FixedVerticesInformation, _T("&Information") );

  wxMenu * menuInterpolate = new wxMenu;
  menuInterpolate->Append( ID_InterpolateLinearModes, _T("&Linear modes...") );
  menuInterpolate->Append( ID_InterpolateNonLinearModes, _T("&Nonlinear modes...") );
  menuInterpolate->Append( ID_InterpolateModalDerivatives, _T("&Modal derivatives...") );
  menuInterpolate->Append( ID_InterpolateSimulationData, _T("&External simulation data...") );
  menuInterpolate->Append( ID_SaveInterpolant, _T("&Export interpolant...") );

  wxMenu * menuMesh = new wxMenu;
  menuMesh->Append( ID_LoadMesh, _T("&Load triangle mesh..."));
  menuMesh->Append( ID_RemoveTriangleMesh, _T("&Remove triangle mesh"));
  menuMesh->Append( ID_Voxelize, _T("&Voxelize...") );
  menuMesh->AppendSubMenu(menuInterpolate, 
    _T("&Interpolate to triangle mesh"), _T("Interpolate various matrices from simulation mesh to triangle mesh"));

  menuMesh->AppendSeparator();
  menuMesh->Append( ID_LoadVoxelization, _T("L&oad simulation mesh..."));
  menuMesh->Append( ID_SaveVoxelization, _T("&Save simulation mesh..."));
  menuMesh->Append( ID_ExportSurfaceMesh, _T("&Export surface of simulation mesh..."));
  menuMesh->Append( ID_MaterialProperties, _T("Set &material properties..."));
    
  menuMesh->AppendSeparator();
  menuMesh->Append( ID_MeshInformation, _T("&Information"));

  wxMenu * menuFrequencies = new wxMenu;
  menuFrequencies->Append(ID_LoadFrequencies, _T("&Load...") );
  menuFrequencies->Append(ID_SaveFrequencies, _T("&Save...") );
  menuFrequencies->Append(ID_DisplayFrequencies, _T("&Display") );
  menuFrequencies->AppendSeparator();
  menuFrequencies->Append(ID_ScaleTo1Hz, _T("&Scale to 1Hz") );
  menuFrequencies->Append(ID_ScaleFrequencies, _T("&Scale...") );
  menuFrequencies->AppendSeparator();
  menuFrequencies->Append(ID_DeleteFrequencies, _T("&Delete...") );

  wxMenu * menuRigidModes = new wxMenu;
  menuRigidModes->Append(ID_DisplayRigidModes, _T("&Display...") );
  menuRigidModes->Append(ID_SetRigidModes, _T("&Set...") );

  wxMenu * menuImportExport = new wxMenu;
  menuImportExport->Append(ID_ConvertTextMatrixToBinaryMatrix, _T("Import matrix from text file..."));
  menuImportExport->Append(ID_ConvertBinaryMatrixToTextMatrix, _T("Export matrix to text file..."));

  wxMenu * menuExportMassMatrix = new wxMenu;
  menuExportMassMatrix->Append( ID_ExportMassMatrixFull, _T("Full...") );
  menuExportMassMatrix->Append( ID_ExportMassMatrixConstrained, _T("With constrained DOFs removed...") );

  wxMenu * menuExportStiffnessMatrix = new wxMenu;
  menuExportStiffnessMatrix->Append( ID_ExportStiffnessMatrixFull, _T("Full...") );
  menuExportStiffnessMatrix->Append( ID_ExportStiffnessMatrixConstrained, _T("With constrained DOFs removed...") );

  wxMenu * menuLinearModes = new wxMenu;
  menuLinearModes->Append( ID_LoadLinearModes, _T("&Load...") );
  menuLinearModes->Append( ID_SaveLinearModes, _T("&Save...") );
  menuLinearModes->Append( ID_ComputeLinearModes, _T("&Compute...") );
  menuLinearModes->AppendSeparator();
  menuLinearModes->AppendSubMenu(menuFrequencies, _T("Frequencies"), _T("Load/Save/Adjust frequencies of the model"));
  menuLinearModes->AppendSeparator();
  menuLinearModes->AppendSubMenu(menuRigidModes, _T("Rigid modes"), _T("View/Set frequencies corresponding to the rigid modes of the model"));
  menuLinearModes->AppendSeparator();
  menuLinearModes->AppendSubMenu(menuExportMassMatrix, _T("Export &mass matrix"), _T("Export &mass matrix"));
  menuLinearModes->AppendSubMenu(menuExportStiffnessMatrix, _T("Export &stiffness matrix"), _T("Export &stiffness matrix"));
  menuLinearModes->AppendSubMenu(menuImportExport, _T("Import/export any matrix"), _T("Import/export any matrix"));

  // modal derivatives menu
  wxMenu * menuModalDerivatives = new wxMenu;
  menuModalDerivatives->Append( ID_LoadModalDerivatives, _T("&Load modal derivatives...") );
  menuModalDerivatives->Append( ID_SaveModalDerivatives, _T("&Save modal derivatives...") );
  menuModalDerivatives->AppendSeparator();
  menuModalDerivatives->Append( ID_ComputeModalDerivatives, _T("&Compute modal derivatives") );

  // sketch menu
  wxMenu * menuSketch = new wxMenu;
  menuSketch->Append( ID_LoadSketchData, _T("&Load deformation data...") );

  // cubic polynomials menu
  wxMenu * menuCubicPolynomials = new wxMenu;
  menuCubicPolynomials->Append( ID_LoadCubicPolynomials, _T("&Load...") );
  menuCubicPolynomials->Append( ID_SaveCubicPolynomials, _T("&Save...") );
  menuCubicPolynomials->AppendSeparator();
  menuCubicPolynomials->Append( ID_ComputeCubicPolynomials, _T("&Compute...") );

  // nonlinear modes menu
  wxMenu * menuNonLinearModes = new wxMenu;
  menuNonLinearModes->Append( ID_LoadNonLinearModes, _T("&Load...") );
  menuNonLinearModes->Append( ID_SaveNonLinearModes, _T("&Save...") );
  menuNonLinearModes->AppendSeparator();
  menuNonLinearModes->AppendSubMenu(menuModalDerivatives, _T("Modal derivatives"), _T("Obtain nonlinear deformation data from modal derivatives"));
  menuNonLinearModes->AppendSubMenu(menuSketch, _T("External simulation"), _T("Obtain nonlinear deformation data from precomputed deformations (such as from an external offline simulator)"));
  menuNonLinearModes->AppendSeparator();
  menuNonLinearModes->Append( ID_ComputeNonLinearModes, _T("&Compute...") );
  menuNonLinearModes->AppendSeparator();
  menuNonLinearModes->AppendSubMenu(menuCubicPolynomials, _T("Cubic polynomials"), _T("Compute reduced force cubic polynomials"));

  // real-time simulation menu
  wxMenu * realTimeSimulation = new wxMenu;
  realTimeSimulation->Append( ID_BatchPreprocessSingleSimulation, _T("&Batch preprocess single model...") );
  realTimeSimulation->Append( ID_BatchPreprocessManySimulations, _T("&Batch preprocess many models...") );
  realTimeSimulation->AppendSeparator();
  realTimeSimulation->Append( ID_PrepareRuntimeSimulation, _T("&Prepare...") );
  realTimeSimulation->Append( ID_TuneRuntimeSimulation, _T("&Tune...") );
  realTimeSimulation->Append( ID_LaunchRuntimeSimulation, _T("&Launch") );

  // help menu
  wxMenu * menuHelp = new wxMenu;
  menuHelp->Append( ID_MouseHelp, _T("&Mouse help...") );
  menuHelp->Append( ID_ExampleHelp, _T("&Example usage...") );
  menuHelp->Append( ID_TuneRuntimeSimulation, _T("&Real-time simulation tuning help...") );
  menuHelp->AppendSeparator();
  menuHelp->Append( ID_About, _T("&About...") );

  menuBar = new wxMenuBar;
  menuBar->Append( menuFile, _T("&File") );
  menuBar->Append( menuView, _T("&View") );
  menuBar->Append( menuMesh, _T("&Mesh") );
  menuBar->Append( menuFixedVertices, _T("F&ixed Vertices") );
  menuBar->Append( menuLinearModes, _T("&Linear Modes") );
  menuBar->Append( menuNonLinearModes, _T("&Simulation Basis") );
  //menuBar->Append( menuCubicPolynomials, _T("&Cubic Polynomials") );
  menuBar->Append( realTimeSimulation, _T("&Real-time Sim") );
  menuBar->Append( menuHelp, _T("&Help") );

  SetMenuBar( menuBar );

  // create toolbar
  toolBar = CreateToolBar(wxNO_BORDER | wxHORIZONTAL, -1);

  viewText = new wxStaticText(toolBar, -1, 
     _T(" View: No model loaded"), wxDefaultPosition, wxDefaultSize, 
     wxALIGN_CENTRE, _T( "staticText"));

  modeSelectionText = new wxStaticText(toolBar, -1, 
     _T("Mode:"), wxDefaultPosition, wxDefaultSize, 
     wxALIGN_CENTRE, _T( "staticText"));

  modeSelectionControl = new wxSpinCtrl(toolBar, ID_ModeSelection, 
    wxEmptyString, wxDefaultPosition, wxSize(50,-1), wxSP_ARROW_KEYS | wxSP_WRAP, 1, 1000000, 1, _T("wxSpinCtrl"));

  frequencyText = new wxStaticText(toolBar, -1, _T("Frequency [cycles/s]:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE, _T( "staticText"));

  frequencyTextCtrl = new wxTextCtrl(toolBar, -1, _T("--"), wxDefaultPosition, wxSize(80,-1), wxTE_READONLY);

  modeAmplitudeText = new wxStaticText(toolBar, -1, _T("Amplitude:"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE, _T( "staticText"));

  modeAmplitudeControl = new wxTextCtrl(toolBar, ID_RenderingAmplitude, _T("--"), wxDefaultPosition, wxSize(50,-1), wxTE_LEFT);

  toolBar->AddControl(viewText);
  toolBar->AddSeparator();
  toolBar->AddControl(modeSelectionText);
  toolBar->AddControl(modeSelectionControl);
  toolBar->AddSeparator();
  toolBar->AddControl(frequencyText);
  toolBar->AddControl(frequencyTextCtrl);
  toolBar->AddSeparator();
  toolBar->AddControl(modeAmplitudeText);
  toolBar->AddControl(modeAmplitudeControl);
  toolBar->Realize();
  EnableModalToolbar(false);

  wxBoxSizer * myFrameSizer = new wxBoxSizer(wxHORIZONTAL);

  // create OpenGL Window
  CreateOpenGLWindow();
  myFrameSizer->Add(myGLCanvas, 2, wxEXPAND);
  this->SetSizer(myFrameSizer);
  CreateStatusBar();

  wxIcon * myIcon = new wxIcon(iconBits);
  SetIcon(*myIcon);

  menuBar->Check(ID_ShowAxes, true);
  menuBar->Check(ID_ShowMaterialGroups, true);

  SelectView(uiState.viewMode);
  modeSelectionControl->SetValue(1);

  UpdateMenus();

  myGLCanvas->Reshape();
  myGLCanvas->Refresh(false);

  uiState.vertexSelectionActivated = false;
  uiState.enableBatchOutput = true;
  uiState.renderMesh = true;

  // init threads
  computationRunning = 0;

  printf(" done.\n");
}

void MyFrame::OnMenuClose(wxMenuEvent & event)
{
  myGLCanvas->Refresh(false);
}

void MyFrame::OnMenuOpen(wxMenuEvent & event)
{
  myGLCanvas->Refresh(false);
}

void MyFrame::OnLoadProject(wxCommandEvent& event)
{
}

void MyFrame::OnSaveProject(wxCommandEvent& event)
{
}

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
  Close(TRUE);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
  wxString message = 
    _T( "Modal Deformation Precomputation Utility.\n"
        L"Version ") + versionString + _T("\n") +
    _T( "\n"
        L"Computes:\n"
        L"- Linear vibrational modes (both constrained and free models are supported)\n"
        L"- Modal derivatives\n"
        L"- Simulation basis (using linear modes and modal derivatives, or external simulation data)\n"
        L"- Reduced cubic polynomials for geometrically nonlinear simulations\n"
        L"- A cubic (voxel) mesh to the input triangle geometry\n"
        L"\n"
        L"By Jernej Barbic, CMU, MIT, USC, 2007-2013.\n"
  );

  wxMessageBox(message,
    _T("About Large Modal Deformation Factory"), wxOK | wxICON_INFORMATION, this);
}

void MyFrame::OnMouseHelp(wxCommandEvent& WXUNUSED(event))
{
  wxString message = 
  _T( "LEFT mouse button:\n"
  L"  click: constrain a vertex\n"
  L"  shift-click: toggle constraining a vertex\n"
  L"  drag: marquee-select constrained vertices\n"
  L"  shift-drag: marquee-toggle constraining vertices\n"
  L"  NOTE: The \"Fixed Vertices.Select\" menu flag must be activated\n"
  L"\n" 
  L"MIDDLE mouse button:\n"
  L"  drag: camera zoom in/out\n"
  L"\n"
  L"RIGHT mouse button:\n"
  L"  drag: change camera viewpoint\n");
        
  wxMessageBox(message,
    _T("Mouse Help"), wxOK | wxICON_INFORMATION, this);
}

void MyFrame::OnExampleHelp(wxCommandEvent& WXUNUSED(event))
{
  wxString message = 
    _T( "Example usage\n"
        L"Note: This text was also printed to your command prompt window. It is also available in \"instructions.txt\" (shortcut is in the Start menu).\n"
        L"=========================================================================================\n"
        L"\n"
        L"1. Load \"examples/simpleBridge.obj\" (using Mesh.Load triangle mesh)\n\n"
        L"2. Create a voxel simulation mesh (using Mesh.Voxelize, keep the resolution at the default (65))\n\n"
        L"3. (optional) Inspect the properties of the newly created simulation mesh (using Mesh.Information)\n\n"
        L"4. (optional) Switch the view to the triangle mesh (and later back to simulation mesh) using View.Triangle Mesh, and View.Simulation Mesh\n\n"
        L"5. (optional) Select some vertices that are to be permanently fixed (left-click with the mouse, use shift to add/subtract, you can also marquee select)\n\n"
        L"6. (optional) At any time, you can inspect the program text output in the text window (might be initially hidden behind the main program window)\n\n"
        L"7. Compute the linear modes (using Linear Modes.Compute; compute, say, 10 modes)\n\n"
        L"8. (optional) See the different modes and their frequencies by pressing the up arrow/down arrow on the toolbar\n\n"
        L"9. Compute modal derivatives (using the Simulation Basis.Modal derivatives.Compute modal derivatives)\n\n"
        L"10. Compute nonlinear modes (using Simulation Basis.Compute; compute, say, 15 modes)\n\n"
        L"11. Generate cubic polynomials (using Cubic polynomials.Compute; this part of the computation is multi-threaded, so it will work faster on a multi-core computer)\n\n"
        L"12. Prepare the real-time simulation (using Real-time Sim.Prepare real-time simulation)\n\n"
        L"13. Launch the real-time simulation (using Real-time Sim.Launch real-time simulation). Note: reducedDynamicSolver-rt must be in the path, check the .bat file.\n\n"
        L"14. Interact with the real-time simulation using the mouse (left-click on the model and pull on it)\n\n"
        L"15. Close the real-time simulation window\n" );
        
  printf("%s", (const char*)message.mb_str());

  wxMessageBox(message, _T("Example usage"), wxOK | wxICON_INFORMATION, this);
}

void MyFrame::CreateOpenGLWindow()
{
  int attribList[32];
  attribList[0] = WX_GL_RGBA;
  attribList[1] = WX_GL_DOUBLEBUFFER;
  attribList[2] = WX_GL_DEPTH_SIZE;
  attribList[3] = 16;
  attribList[4] = WX_GL_STENCIL_SIZE;
  attribList[5] = 8;
  attribList[6] = 0;

  myGLCanvas = new MyGLCanvas(&precomputationState, &uiState, 
    this, -1, 
    wxPoint(0,0), wxSize(640,640),
    wxSUNKEN_BORDER, _T("Visualization"), attribList);
  myGLCanvas->Refresh(false);
}

void MyFrame::InitOpenGL()
{
  myGLCanvas->Show();
  //myGLCanvas->SetCurrent();
  myGLCanvas->InitOpenGL();
  wxSafeYield();
}

void MyFrame::StopComputation(int parameter)
{
}

void MyFrame::SpawnComputationThread(void * (MyFrame::*ComputationThread) (void*))
{
}

void MyFrame::errMsg(wxString title, wxString message, wxWindow* parent)
{
  WxDialogWrapper<wxMessageDialog> dlg( new wxMessageDialog(parent,
    wxString((const char*)message.mb_str(), wxConvUTF8), 
    wxString((const char*)title.mb_str(), wxConvUTF8),
    wxOK | wxICON_ERROR) );
  dlg.get()->ShowModal();
}

void MyFrame::OnMenu(wxCommandEvent& event)
{
  myGLCanvas->Refresh(FALSE);
}

