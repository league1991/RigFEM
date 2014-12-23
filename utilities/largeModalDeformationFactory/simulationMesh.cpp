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

// manage the simulation mesh

#include "volumetricMeshLoader.h"
#include "generateSurfaceMesh.h"
#include "volumetricMeshENuMaterial.h"
#include "objMeshOffsetVoxels.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnVoxelize(wxCommandEvent& event)
{
  char cubicMeshResolutionStringC[256];
  sprintf(cubicMeshResolutionStringC, "%d", uiState.cubicMeshResolution);

  string cubicMeshResolutionString(cubicMeshResolutionStringC);

  wxDialog * dlg = new wxDialog(this, -1, _T("Select resolution"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * resolutionText = new wxStaticText(dlg, -1, 
       _T("Voxelization resolution: "), wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER, _T( "staticText"));

  int maxResolution = 16384;
  wxSpinCtrl * resolutionControl = new wxSpinCtrl(dlg, -1, 
      wxEmptyString, wxDefaultPosition, wxSize(70,-1), wxSP_ARROW_KEYS, 
      1, maxResolution, uiState.cubicMeshResolution, _T("wxSpinCtrl"));
  resolutionControl->SetValue(uiState.cubicMeshResolution);

  wxBoxSizer * resolutionSelectorSizer = new wxBoxSizer(wxHORIZONTAL);
  resolutionSelectorSizer->Add(resolutionText, 0, wxCENTER);
  resolutionSelectorSizer->Add(resolutionControl, 0, wxCENTER);
  dlgSizer->Add(resolutionSelectorSizer, 0, wxALL | wxCENTER, 8);

  wxCheckBox * interiorChamberCheckBox = new wxCheckBox(dlg, 
    -1, _T("Fill interior chambers"));
  interiorChamberCheckBox->SetValue(uiState.fillInteriorChambers);
  dlgSizer->Add(interiorChamberCheckBox, 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  int value = resolutionControl->GetValue();
  uiState.fillInteriorChambers = interiorChamberCheckBox->GetValue();

  delete(dlg);

  bool goodInput = true;
  if (goodInput)
  {
    uiState.cubicMeshResolution = value;

    SetCursor(*wxHOURGLASS_CURSOR);
    VoxelizationWorker(precomputationState.renderingMesh, 
      uiState.cubicMeshResolution, uiState.fillInteriorChambers);
    SetCursor(*wxSTANDARD_CURSOR);
printf("Worker done.\n");fflush(NULL);

    precomputationState.simulationMeshAvailable = true;
    precomputationState.fixedVertices.clear();
    ActivateVertexSelection(true);

    ComputeSimulationMeshGeometricParameters();
    SetDefaultMaterialParameters();

    DeallocateSimulationData();

    UpdateMenus();

    myGLCanvas->UpdateSimulationMeshRenderData();
    SelectView(UIState::VIEW_SIMULATION_MESH);
  }
  else
  {
    this->errMsg( _T("Invalid resolution"),  _T("Invalid resolution: "));
  }
}

int MyFrame::LoadSimulationMesh(wxString & meshFilename)
{
  VolumetricMesh * mesh;
  try
  {
    SetCursor(*wxHOURGLASS_CURSOR);
    const char * filename = meshFilename.mb_str();
    mesh = VolumetricMeshLoader::load((char*)filename);
    SetCursor(*wxSTANDARD_CURSOR);
  }
  catch( int exceptionCode )
  {
    exceptionCode++; // to avoid the warning
    this->errMsg( _T("Loading error"),
      _T("Unable to load simulation mesh from ") +  meshFilename );
    return 1;
  }

  // check if the mesh consists of ENu materials
  bool ENuOnly = true;
  for(int i=0; i<mesh->getNumMaterials(); i++)
  {
    VolumetricMesh::Material * material = mesh->getMaterial(i);
    VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
    if (eNuMaterial == NULL)
      ENuOnly = false;
  }

  if (!ENuOnly)
  {
    this->errMsg( _T("Loading error"), _T("Simulation mesh in ") +  meshFilename + _T(" does not consist only of ENU materials."));
    return 2;
  }

  delete(precomputationState.simulationMesh);
  precomputationState.simulationMesh = mesh;
  uiState.projectName = meshFilename.AfterLast('/').AfterLast('\\').BeforeLast('.');

  precomputationState.simulationMeshAvailable = true;
  precomputationState.fixedVertices.clear();
  ActivateVertexSelection(true);

  ComputeSimulationMeshGeometricParameters();

  DeallocateSimulationData();

  UpdateMenus();

  myGLCanvas->UpdateSimulationMeshRenderData();
  SelectView(UIState::VIEW_SIMULATION_MESH);

  return 0;
}

void MyFrame::OnLoadSimulationMesh(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load simulation mesh (cubic or tetrahedral)"), uiState.currentWorkingDirectory, _T(""), _T("Vega Mesh Files(*.veg)|*.veg|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString meshFilename = dlg->GetPath();
    dlg->Destroy();
    SaveCurrentWorkingDirectory(meshFilename);
  
    if( !meshFilename.empty() )
    {
      int code = LoadSimulationMesh(meshFilename);
      if (code != 0)
      {
        printf("Error loading the simulation mesh from %s, exit code %d.\n", (const char*)meshFilename.mb_str(), code);
        return;
      }
    }
  }
  else
  {
    dlg->Destroy();
  }
}

void MyFrame::OnSaveSimulationMesh(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Save cubic simulation mesh"), uiState.currentWorkingDirectory, _T(""), _T("Vega Mesh Files(*.veg)|*.veg|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString meshFilename = dlg->GetPath();
    SaveCurrentWorkingDirectory(meshFilename);
    if( !meshFilename.empty() )
    {
      const char * filename = meshFilename.mb_str();
      int code = precomputationState.simulationMesh->saveToAscii((char*)filename);
      if (code != 0)
        this->errMsg( _T("Saving error"), 
          _T("Unable to save simulation mesh to ") +  meshFilename );
    }
  }

  dlg->Destroy();
}

void MyFrame::OnMeshInformation(wxCommandEvent& event)
{
  string title = string("Mesh Information");

  TCHAR meshString[4096];

  if (precomputationState.renderingMeshAvailable)
  {
    int numVertices = (int)precomputationState.renderingMesh->getNumVertices();
    int numFaces = (int)precomputationState.renderingMesh->getNumFaces();
    int numGroups = (int)precomputationState.renderingMesh->getNumGroups();
  
    wsprintf(meshString, 
      L"  Num vertices: %d  \n"
      L"  Num faces: %d  \n"
      L"  Num groups: %d  \n",
      numVertices, numFaces, numGroups);  
  }
  else
  {
    wsprintf(meshString, L"  Not available.\n");
  }

  TCHAR simulationString[4096];

  if (precomputationState.simulationMeshAvailable)
  {
    wsprintf(simulationString, 
      L"  Num vertices: %d  \n"
      L"  Num elements: %d  \n"
      L"  Element type: %s  \n"
      L"  Num materials: %d  \n"
      L"  Num regions: %d  \n",
      precomputationState.simulationMesh->getNumVertices(),
      precomputationState.simulationMesh->getNumElements(),
      (precomputationState.simulationMesh->getElementType() == VolumetricMesh::TET) ? "TET" : (precomputationState.simulationMesh->getElementType() == VolumetricMesh::CUBIC ? "CUBIC" : "UNKNOWN"),
      precomputationState.simulationMesh->getNumMaterials(),
      precomputationState.simulationMesh->getNumRegions()
      );
  }
  else
  {
    wsprintf(simulationString, L"  Not available.  \n");
  }

  char information[9192];
  sprintf(information, "Triangle mesh:\n%s\nSimulation mesh:\n%s", meshString, simulationString);

  wxWindow * parent = NULL;
  WxDialogWrapper<wxMessageDialog> dlg( new wxMessageDialog(parent,
    wxString((const char*)information, wxConvUTF8), 
    wxString(title.c_str(), wxConvUTF8), wxOK | wxICON_INFORMATION) );
  dlg.get()->ShowModal();
}

void * MyFrame::VoxelizationWorker(ObjMesh * objMesh, int resolution, bool fillHoles)
{
  double expansionFactor = 1.2;
  int depth = 0;
  int resolutionArray[3] = { resolution, resolution, resolution };
  ObjMeshOffsetVoxels * objMeshOffsetVoxels = 
    new ObjMeshOffsetVoxels( objMesh, resolutionArray, depth, expansionFactor );

  if (fillHoles)
  {
     vector<Vec3d> componentSeeds; 
     vector<int> componentSize;
     objMeshOffsetVoxels->emptyComponents(componentSeeds, componentSize);
     objMeshOffsetVoxels->floodFill(componentSeeds);
  }

  int numVertices;
  int numElements;
  double * vertices;
  int * voxels;

  free(precomputationState.simulationSurfaceMesh);
  free(precomputationState.FFDinterpolationVertices);
  free(precomputationState.FFDinterpolationWeights);

  objMeshOffsetVoxels->generateCubicMesh(
      &numVertices, &vertices, 
      &numElements, &voxels,
      &precomputationState.FFDinterpolationVertices, 
      &precomputationState.FFDinterpolationWeights,  
      &precomputationState.simulationSurfaceMesh);

  delete(objMeshOffsetVoxels);


  double E = 1E9;
  double nu = 0.45;
  double density = 1000.0;
  delete(precomputationState.simulationMesh);
  precomputationState.simulationMesh = new CubicMesh(numVertices, vertices, numElements, voxels, E, nu, density);

  printf("Generated the cubic simulation mesh. Vertices: %d | Voxels: %d\n", numVertices, numElements);

  free(vertices);
  free(voxels);

  computationRunning = -1;

  return NULL;
}

void * MyFrame::VoxelizationThread(void * vargp)
{
  return NULL;
}

void MyFrame::ComputeSimulationMeshGeometricParameters()
{
  Vec3d center;
  precomputationState.simulationMesh->getMeshGeometricParameters
    (center, &(precomputationState.simulationMeshGeometricParameters.radius));

  precomputationState.simulationMeshGeometricParameters.centerX = center[0];
  precomputationState.simulationMeshGeometricParameters.centerY = center[1];
  precomputationState.simulationMeshGeometricParameters.centerZ = center[2];
}

void MyFrame::OnMaterialProperties(wxCommandEvent & event)
{
  if (precomputationState.simulationMesh->getNumMaterials() > 1)
  {
    wxMessageDialog * confirmationDialog = new wxMessageDialog
      (this, _T(
        "Warning: existing simulation mesh has more than one material. "
        L"This dialog can only edit one material; other materials will be deleted. "
        L"Do you want to continue?"
        ),
      _T("More than one material encountered"), 
      wxYES_NO | wxICON_QUESTION);

    if (confirmationDialog->ShowModal() != wxID_YES)
    {
      delete(confirmationDialog);
      return;
    }

    delete(confirmationDialog);
  }

  wxDialog * dlg = new wxDialog(this, -1, _T("Material properties"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  // Young's modulus
  char EString[96];
  if (precomputationState.simulationMesh->getNumMaterials() >= 1)
  {
    VolumetricMesh::Material * material = precomputationState.simulationMesh->getMaterial(0);
    VolumetricMesh::ENuMaterial * eNuMaterial =  downcastENuMaterial(material);
    double value = eNuMaterial->getE();

    int magnitudeOrder = (int)(log(value) / log(10.0));
    int numDigits;
    if (magnitudeOrder > 6)
      numDigits = 12 - magnitudeOrder;
    else
      numDigits = 6;
    if (numDigits < 0)
      numDigits = 0;
    char formatString[16];
    sprintf(formatString, "%%.%df", numDigits);
    sprintf(EString, formatString, value);
  }
  else
    EString[0] = 0;

  wxStaticText * EText = new wxStaticText(dlg, -1, 
       _T("Young's modulus [N/m2]:"), wxDefaultPosition, wxDefaultSize, 
       wxTE_LEFT, _T("staticText"));

  wxTextCtrl * EControl = new wxTextCtrl(dlg, -1, 
      wxString(EString, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * ESizer = new wxBoxSizer(wxHORIZONTAL);
  ESizer->Add(EText, 0, wxCENTER);
  ESizer->Add(EControl, 0, wxCENTER);
  dlgSizer->Add(ESizer, 0, wxLEFT | wxTOP | wxRIGHT, 8);

  // Poisson's ratio
  char nuString[96];
  if (precomputationState.simulationMesh->getNumMaterials() >= 1)
  {
    VolumetricMesh::Material * material = precomputationState.simulationMesh->getMaterial(0);
    VolumetricMesh::ENuMaterial * eNuMaterial =  downcastENuMaterial(material);
    double value = eNuMaterial->getNu();
    sprintf(nuString, "%f", value);
  }
  else
    nuString[0] = 0;

  wxStaticText * nuText = new wxStaticText(dlg, -1, 
       _T("Poisson ratio:"), wxDefaultPosition, wxDefaultSize, 
       wxTE_LEFT, _T( "staticText"));

  wxTextCtrl * nuControl = new wxTextCtrl(dlg, -1, 
      wxString(nuString, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * nuSizer = new wxBoxSizer(wxHORIZONTAL);
  nuSizer->Add(nuText, 0, wxCENTER);
  nuSizer->Add(nuControl, 0, wxCENTER);
  dlgSizer->Add(nuSizer, 0, wxLEFT | wxTOP | wxRIGHT, 8);

  // density
  char rhoString[96];
  if (precomputationState.simulationMesh->getNumMaterials() >= 1)
  {
    double value = precomputationState.simulationMesh->getMaterial(0)->getDensity();
    sprintf(rhoString, "%f", value);
  }
  else
    rhoString[0] = 0;

  wxStaticText * rhoText = new wxStaticText(dlg, -1, 
       _T("Mass density [kg/m3]:"), wxDefaultPosition, wxDefaultSize, 
       wxTE_LEFT, _T( "staticText"));

  wxTextCtrl * rhoControl = new wxTextCtrl(dlg, -1, 
      wxString(rhoString, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * rhoSizer = new wxBoxSizer(wxHORIZONTAL);
  rhoSizer->Add(rhoText, 0, wxCENTER);
  rhoSizer->Add(rhoControl, 0, wxCENTER);
  dlgSizer->Add(rhoSizer, 0, wxALL, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 
    0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  double E, nu, density;
  
  bool EError = !(EControl->GetValue().ToDouble(&E));
  EError = EError || (E <= 0);
  if (EError)
  {
    errMsg( _T("Invalid Young's modulus"), _T("Invalid Young's modulus: ") + EControl->GetValue());
    return;
  }

  bool nuError = !(nuControl->GetValue().ToDouble(&nu));
  nuError = nuError || (nu <= -0.5) || (nu >= 0.5);
  if (nuError)
  {
    errMsg(_T("Invalid Poisson ratio"),
      _T("Invalid Poisson ratio: ") + nuControl->GetValue());
    return;
  }

  bool rhoError = !(rhoControl->GetValue().ToDouble(&density));
  rhoError = rhoError || (density <= 0);
  if (rhoError)
  {
    errMsg(_T("Invalid mass density"),
      _T("Invalid mass density: ") + rhoControl->GetValue());
    return;
  }

  printf("Setting: E=%G nu=%G density=%G\n", E, nu, density);
  precomputationState.simulationMesh->setSingleMaterial(E, nu, density);

  delete(dlg);   
}

void MyFrame::OnExportSurfaceMesh(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Export surface mesh of the simulation mesh"), uiState.currentWorkingDirectory, _T(""), _T("Obj Files(*.obj)|*.obj|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString objMeshname(dlg->GetPath().GetData());
    SaveCurrentWorkingDirectory(objMeshname);
    if( !objMeshname.empty() )
    {
      const char * filename = objMeshname.mb_str();
      ObjMesh * newObjMesh;
      SetCursor(*wxHOURGLASS_CURSOR);
      try
      {
        newObjMesh = GenerateSurfaceMesh::ComputeMesh(precomputationState.simulationMesh);
        newObjMesh->save(string(filename));
        delete(newObjMesh);
      }      
      catch( ObjMeshException myException )
      {
        SetCursor(*wxSTANDARD_CURSOR);
        this->errMsg( _T("Error creating the surface mesh"),  wxString(myException.getReason().c_str(), wxConvUTF8) );
        dlg->Destroy();
        return;
      }
      SetCursor(*wxSTANDARD_CURSOR);
    }
  }
  dlg->Destroy();
}

void MyFrame::SetDefaultMaterialParameters()
{
  double E, nu, density;
  precomputationState.simulationMesh->getDefaultMaterial(&E, &nu, &density);
  precomputationState.simulationMesh->setSingleMaterial(E, nu, density);
}

void MyFrame::DeleteInterpolationData()
{
  precomputationState.interpolationDataAvailable = false;
  free(precomputationState.interpolationData_vertices);
  precomputationState.interpolationData_vertices = NULL;
  free(precomputationState.interpolationData_weights);
  precomputationState.interpolationData_weights = NULL;
}

void MyFrame::DeallocateSimulationData()
{
  // fixed vertices
  precomputationState.fixedVerticesAvailable = false;
  precomputationState.fixedVertices.clear();
  
  // linear modes
  precomputationState.linearModesAvailable = false;
  delete(precomputationState.linearModalMatrix);
  precomputationState.linearModalMatrix = NULL;
  precomputationState.rLin = 0;

  // num rigid modes
  precomputationState.numRigidModes = 0;

  // frequencies
  precomputationState.frequenciesAvailable = false;
  free(precomputationState.frequencies);
  precomputationState.frequencies = NULL;

  // nonlinear modes
  precomputationState.nonLinearModesAvailable = false;
  delete(precomputationState.nonLinearModalMatrix);
  precomputationState.nonLinearModalMatrix = NULL;
  precomputationState.rNonLin = 0;

  // modal derivatives
  precomputationState.modalDerivativesAvailable = false;
  delete(precomputationState.modalDerivativesMatrix);
  precomputationState.modalDerivativesMatrix = NULL;

  // sketch
  precomputationState.sketchDataAvailable = false;
  delete(precomputationState.sketchDataMatrix);
  precomputationState.sketchDataMatrix = NULL;

  // cubic polynomials
  precomputationState.cubicPolynomialsAvailable = false;
  delete(precomputationState.cubicPolynomials);
  precomputationState.cubicPolynomials = NULL;

  DeleteInterpolationData();

  precomputationState.runtimeSimReady = false;

  UpdateMenus();
}

