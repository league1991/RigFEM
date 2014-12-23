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

// linear mode computation

#include <math.h>
#include <float.h>
#include "sparseMatrix.h"
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "StVKCubeABCD.h"
#include "matrixIO.h"
#include "insertRows.h"
#include "StVKElementABCDLoader.h"
#include "ARPACKSolver.h"
#include "largeModalDeformationFactory.h"

#ifdef WIN32
  #define ISNAN _isnan
#else
  #define ISNAN isnan
#endif

void MyFrame::SetAutoRenderingMagnitude(ModalMatrix * modalMatrix)
{
  // determine rendered magnitude
  // renderedMagnitude * maxAmplitude = 0.1 * radius
  double maxMagnitude = MaxModeMagnitude(modalMatrix);
  double renderingMagnitude;

  if (maxMagnitude < -0.1)
  {
    renderingMagnitude = 1.0;    
          
    wxMessageDialog * dlg = new wxMessageDialog
          (this, _T("Warning: the modes contain NaN entries."),
          _T("NaN entries detected in the modes"), 
           wxOK | wxICON_EXCLAMATION);

    dlg->ShowModal();
    dlg->Destroy();
  }
  else if (maxMagnitude < 1e-8)
  {
    renderingMagnitude = 1.0;    
          
    wxMessageDialog * dlg = new wxMessageDialog (this, _T("Warning: the modes have extremely small magnitude (might be zero)."), _T("Extremely small magnitude modes"), wxOK | wxICON_EXCLAMATION);

    dlg->ShowModal();
    dlg->Destroy();
  }
  else 
    renderingMagnitude = 0.35 * precomputationState.simulationMeshGeometricParameters.radius / maxMagnitude;

  char s[4096];
  sprintf(s, "%.3f", renderingMagnitude);
  modeAmplitudeControl->SetValue(wxString(s, wxConvUTF8));

  UpdateModalToolbar();
}

double MyFrame::MaxModeMagnitude(ModalMatrix * modalMatrix)
{
  double * U = modalMatrix->GetMatrix();
  int n = modalMatrix->Getn();
  int r = modalMatrix->Getr();

  double maxMagnitude = 0;
  for(int i=0; i< n*r; i++)
  {
    double * entry = &U[3*i];
    if (ISNAN(entry[0]) || ISNAN(entry[1]) || ISNAN(entry[2]))
      return -1;
    double magnitude = sqrt(entry[0]*entry[0] + entry[1]*entry[1] + entry[2]*entry[2]);
    if (magnitude > maxMagnitude)
      maxMagnitude = magnitude;
  }

  return maxMagnitude;
}

void MyFrame::OnLoadLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load linear modes"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.Ulin)|*.Ulin|Modal Matrix Files(*.U)|*.U|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString linearModesFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(linearModesFilename);
    if( !linearModesFilename.empty() )
    {
      int newr;
      double * newLinearModes = NULL;

      int n1;
      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = linearModesFilename.mb_str();
      int code = ReadMatrixFromDisk((char*)filename, &n1, &newr, &newLinearModes);
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  _T("Unable to load linear modes from ") + linearModesFilename );
        dlg->Destroy();
        return;
      }

      if (n1 != 3 * precomputationState.simulationMesh->getNumVertices())
      {
        this->errMsg( _T("Loading error"), _T("The number of vertices in ") + linearModesFilename + _T(" does not match the simulation mesh."));
        free(newLinearModes);
        dlg->Destroy();
        return;
      }

      if (precomputationState.frequenciesAvailable)
      {
        // check that the number of modes is consistent with the existing number of frequencies
        if (newr != precomputationState.rLin)
        {
          wxMessageDialog * confirmationDialog = new wxMessageDialog (this, _T("Warning: number of existing frequencies does not match the number of modes. Delete existing frequencies?"), _T("Mismatch in the number of frequencies"), wxYES_NO | wxICON_EXCLAMATION);

          if (confirmationDialog->ShowModal() != wxID_YES)
          {
            free(newLinearModes);
            delete(confirmationDialog);
            dlg->Destroy();
            return;
          }
          else
          {
            delete(confirmationDialog);
            free(precomputationState.frequencies);
            precomputationState.frequenciesAvailable = false;
          }
        }
      }

      // success
      delete(precomputationState.linearModalMatrix);
      precomputationState.rLin = newr;
      precomputationState.linearModalMatrix = new ModalMatrix(
        precomputationState.simulationMesh->getNumVertices(), precomputationState.rLin, newLinearModes);
      free(newLinearModes);

      precomputationState.linearModesAvailable = true;

      uiState.numComputedNonLinearModes = 2 * (precomputationState.rLin - precomputationState.numRigidModes);
      uiState.eraseRangeHi = precomputationState.rLin;

      modeSelectionControl->SetValue(1);

      SelectView(UIState::VIEW_LINEAR_MODES);
      SetAutoRenderingMagnitude(precomputationState.linearModalMatrix);

      UpdateMenus();

      myGLCanvas->UpdateLinearModesRenderData();

      Refresh();
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSaveLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Save linear modes"),
    uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.Ulin)|*.Ulin|All files(*.*)|*.*"),
    wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString linearModesFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(linearModesFilename);
    if( !linearModesFilename.empty() )
    {
      const char * filename = linearModesFilename.mb_str();
      int code = WriteMatrixToDisk((char*)filename,
        3 * precomputationState.linearModalMatrix->Getn(), 
        precomputationState.linearModalMatrix->Getr(), 
        precomputationState.linearModalMatrix->GetMatrix());

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save linear modes to ") + linearModesFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnExportLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Export linear modes"),
	uiState.currentWorkingDirectory, _T(""), _T("Text Files(*.txt)|*.txt|All files(*.*)|*.*"),
	wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);
  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString linearModesFilename( dlg->GetPath().GetData() );
    SaveCurrentWorkingDirectory(linearModesFilename);
    if( !linearModesFilename.empty() )
    {
      const char * filename = linearModesFilename.mb_str();
      FILE * fout = fopen((char*)filename, "w");
      if (fout)
      {
        double * U = precomputationState.linearModalMatrix->GetMatrix();
        int n = precomputationState.linearModalMatrix->Getn();
        fprintf(fout,"%d\n%d\n", 3*n, precomputationState.rLin);
        for(int i=0; i<3*n; i++)
        {
          for(int j=0; j<precomputationState.rLin; j++)
          {
            fprintf(fout, "%.15f ", U[ELT(3*n,i,j)]);
          }
          fprintf(fout,"\n");
        }
        fclose(fout);
      }
      else
      {
        this->errMsg( _T("Exporting error"),  
          _T("Unable to export linear modes to ") + linearModesFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnComputeLinearModes(wxCommandEvent& event)
{
  char numLinearModesStringC[256];

  int numSuggestedModes = uiState.numComputedLinearModes;
  if (uiState.firstModalComputation)
  {
    // increase the number of suggested modes based on the number of constrained vertices
    if (precomputationState.fixedVertices.size() == 2)
    {
      numSuggestedModes += 1;
    }
    else if (precomputationState.fixedVertices.size() == 1)
    {
      numSuggestedModes += 3;
    }
    else if (precomputationState.fixedVertices.size() == 0)
    {
      numSuggestedModes += 6;
    }
  }

  sprintf(numLinearModesStringC, "%d", numSuggestedModes);

  wxDialog * dlg = new wxDialog(this, -1, _T("Select number of modes"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * numModesText = new wxStaticText(dlg, -1, 
       _T("Number of modes: "), wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numModesControl = new wxTextCtrl(dlg, -1, 
      wxString(numLinearModesStringC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numModesSizer = new wxBoxSizer(wxHORIZONTAL);
  numModesSizer->Add(numModesText, 0, wxCENTER);
  numModesSizer->Add(numModesControl, 0, wxCENTER);
  dlgSizer->Add(numModesSizer, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString valueString = numModesControl->GetValue();

  delete(dlg);

  long value;
  bool goodInput = valueString.ToLong(&value);

  if (goodInput)
  {
    if ((value <= 0) || (value > 16384))
      goodInput = false;
  }

  if (goodInput)
  {
    int oldNumModes = uiState.numComputedLinearModes;
    uiState.numComputedLinearModes = value;

    double * newFrequencies = NULL;
    double * newLinearModes = NULL;
    int newr;

    SetCursor(*wxHOURGLASS_CURSOR);
    LinearModesWorker(
      uiState.numComputedLinearModes,
      &newr, &newFrequencies, &newLinearModes );
    SetCursor(*wxSTANDARD_CURSOR);

    if (newr != uiState.numComputedLinearModes)
    {
      uiState.numComputedLinearModes = oldNumModes;
      this->errMsg( _T("Linear mode computation failed"), 
             _T("Linear mode computation failed.") );
      free(newFrequencies);
      free(newLinearModes);
      return;
    }
    else
    {
      precomputationState.linearModesAvailable = true;
      uiState.firstModalComputation = false;
      UpdateMenus();
 
      delete(precomputationState.linearModalMatrix);
      precomputationState.rLin = uiState.numComputedLinearModes;
      precomputationState.linearModalMatrix = new ModalMatrix(
        precomputationState.simulationMesh->getNumVertices(), 
        precomputationState.rLin, newLinearModes);
      free(newLinearModes);

      free(precomputationState.frequencies);
      precomputationState.frequencies = newFrequencies;
      precomputationState.frequenciesAvailable = true;

      precomputationState.linearModesAvailable = true;
      uiState.eraseRangeHi = precomputationState.rLin;

      // set the rigid modes based on the number of constrained vertices
      if (precomputationState.fixedVertices.size() >= 3)
      {
        precomputationState.numRigidModes = 0;
      }
      else if (precomputationState.fixedVertices.size() == 2)
      {
        precomputationState.numRigidModes = 1;
      }
      else if (precomputationState.fixedVertices.size() == 1)
      {
        precomputationState.numRigidModes = 3;
      }
      else
      {
        precomputationState.numRigidModes = 6;
      }

      uiState.numComputedNonLinearModes = 2 * (precomputationState.rLin - precomputationState.numRigidModes);

      modeSelectionControl->SetValue(1);

      SelectView(UIState::VIEW_LINEAR_MODES);
      SetAutoRenderingMagnitude(precomputationState.linearModalMatrix);

      UpdateMenus();

      myGLCanvas->UpdateLinearModesRenderData();

      Refresh();
    }
  }
  else
  {
    this->errMsg( _T("Invalid number of linear modes"),  _T("Invalid number of linear modes: ") +  valueString );
  }
}

void * MyFrame::LinearModesWorker(
      int numDesiredModes,
      int * r, double ** frequencies_, double ** modes_ )
{
  *r = -1;

  // create mass matrix
  SparseMatrix * massMatrix;
  GenerateMassMatrix::computeMassMatrix(precomputationState.simulationMesh, &massMatrix, true);

  // create stiffness matrix
  StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(precomputationState.simulationMesh);
  StVKInternalForces * internalForces = 
    new StVKInternalForces(precomputationState.simulationMesh, precomputedIntegrals);

  SparseMatrix * stiffnessMatrix;
  StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
  stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
  double * zero = (double*) calloc(3 * precomputationState.simulationMesh->getNumVertices(), sizeof(double));
  stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);

  free(zero);
  delete(precomputedIntegrals);
  delete(stiffnessMatrixClass);
  delete(internalForces);

  // constrain the degrees of freedom
  int numConstrainedVertices = (int) (precomputationState.fixedVertices.size());
  int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
  set<int> :: iterator iter;
  int i = 0;
  for(iter = precomputationState.fixedVertices.begin(); iter != precomputationState.fixedVertices.end(); iter++)
  {
    constrainedDOFs[3*i+0] = 3 * (*iter) + 1;
    constrainedDOFs[3*i+1] = 3 * (*iter) + 2;
    constrainedDOFs[3*i+2] = 3 * (*iter) + 3;
    i++;
  }

  int oneIndexed = 1;
  massMatrix->RemoveRowsColumns(
    3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  stiffnessMatrix->RemoveRowsColumns(
    3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  // call ARPACK

  double * frequenciesTemp = (double*) malloc (sizeof(double) * numDesiredModes);
  int numRetainedDOFs = stiffnessMatrix->Getn();
  double * modesTemp = (double*) malloc 
    (sizeof(double) * numDesiredModes * numRetainedDOFs);

  printf("Computing linear modes using ARPACK: ...\n");
  PerformanceCounter ARPACKCounter;
  double sigma = -1.0;

  int numLinearSolverThreads = wxThread::GetCPUCount();
  if (numLinearSolverThreads > 3)
    numLinearSolverThreads = 3; // diminished returns in solver beyond 3 threads

  //massMatrix->Save("MFactory");
  //stiffnessMatrix->Save("KFactory");

  ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv
    (stiffnessMatrix, massMatrix, 
     numDesiredModes, frequenciesTemp, 
     modesTemp, sigma, numLinearSolverThreads);

  ARPACKCounter.StopCounter();
  double ARPACKTime = ARPACKCounter.GetElapsedTime();
  printf("ARPACK time: %G s.\n", ARPACKTime); fflush(NULL);

  if (nconv < numDesiredModes)
  {
    free(modesTemp);
    free(frequenciesTemp);
    *r = -3;
    free(constrainedDOFs);
    delete(massMatrix);
    delete(stiffnessMatrix);
    return NULL;
  }

  int n3 = 3 * precomputationState.simulationMesh->getNumVertices();
  *frequencies_ = (double*) calloc (numDesiredModes, sizeof(double));
  *modes_ = (double*) calloc (numDesiredModes * n3, sizeof(double));

  for(int i=0; i<numDesiredModes; i++)
  {
    // insert zero rows into the computed modes
    int oneIndexed = 1;
    InsertRows(n3, &modesTemp[numRetainedDOFs*i], &((*modes_)[n3*i]), 
      3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
  }

  for(int i=0; i<numDesiredModes; i++)
  {
    if (frequenciesTemp[i] <= 0)
      (*frequencies_)[i] = 0.0;
    else
      (*frequencies_)[i] = sqrt((frequenciesTemp)[i]) / (2 * M_PI);
  }
 
  free(modesTemp);
  free(frequenciesTemp);
  free(constrainedDOFs);

  delete(massMatrix);
  delete(stiffnessMatrix);

  *r = numDesiredModes;

  return NULL;
}

void MyFrame::OnChangeRenderedMode(wxSpinEvent& event)
{
  UpdateModalToolbar();
}

void MyFrame::OnChangeRenderingAmplitude(wxCommandEvent& event)
{
  if ((modeAmplitudeControl != NULL) && (modeAmplitudeControl->GetValue() != _T("--")))
    UpdateModalToolbar();
}

void MyFrame::OnExportStiffnessMatrixFull(wxCommandEvent& event)
{
  ExportStiffnessMatrix(true);
}

void MyFrame::OnExportStiffnessMatrixConstrained(wxCommandEvent& event)
{
  ExportStiffnessMatrix(false);
}

void MyFrame::ExportStiffnessMatrix(bool fullMatrix)
{
  if ((!fullMatrix) && (!precomputationState.fixedVerticesAvailable))
  {
    this->errMsg( _T("Error"),  
      _T("No fixed vertices have been specified.") );
    return;
  }

  wxFileDialog *dlg = new wxFileDialog(this, _T("Export stiffness matrix"), uiState.currentWorkingDirectory, _T(""), _T("Stiffness Matrix Files(*.K)|*.K|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);
  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString stiffnessMatrixFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(stiffnessMatrixFilename);
    if( !stiffnessMatrixFilename.empty() )
    {
      // create stiffness matrix
      StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(precomputationState.simulationMesh);
      StVKInternalForces * internalForces = 
        new StVKInternalForces(precomputationState.simulationMesh, precomputedIntegrals);

      SparseMatrix * stiffnessMatrix;
      StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
      stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
      double * zero = (double*) calloc(3 * precomputationState.simulationMesh->getNumVertices(), sizeof(double));
      stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);

      free(zero);
      delete(precomputedIntegrals);
      delete(stiffnessMatrixClass);
      delete(internalForces);

      if (!fullMatrix)
      {
        // constrain the degrees of freedom
        int numConstrainedVertices = (int) (precomputationState.fixedVertices.size());
        int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
        int i = 0;
        for(set<int> :: iterator iter = precomputationState.fixedVertices.begin(); 
          iter != precomputationState.fixedVertices.end(); iter++)
        {
          constrainedDOFs[3*i+0] = 3 * (*iter) + 1;
          constrainedDOFs[3*i+1] = 3 * (*iter) + 2;
          constrainedDOFs[3*i+2] = 3 * (*iter) + 3;
          i++;
        }

        int oneIndexed = 1;
        stiffnessMatrix->RemoveRowsColumns(
          3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
        free(constrainedDOFs);
      }

      const char * filename = stiffnessMatrixFilename.mb_str();
      int code = stiffnessMatrix->Save((char*)filename);

      delete(stiffnessMatrix);

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save stiffness matrix to ") + stiffnessMatrixFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnExportMassMatrixFull(wxCommandEvent& event)
{
  ExportMassMatrix(true);
}

void MyFrame::OnExportMassMatrixConstrained(wxCommandEvent& event)
{
  ExportMassMatrix(false);
}

void MyFrame::ExportMassMatrix(bool fullMatrix)
{
  if ((!fullMatrix) && (!precomputationState.fixedVerticesAvailable))
  {
    this->errMsg( _T("Error"),  
      _T("No fixed vertices have been specified.") );
    return;
  }

  wxFileDialog *dlg = new wxFileDialog(this, _T("Export mass matrix"), uiState.currentWorkingDirectory, _T(""), _T("Mass Matrix Files(*.M)|*.M|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);
  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString massMatrixFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(massMatrixFilename);
    if( !massMatrixFilename.empty() )
    {
      // create mass matrix

      SparseMatrix * massMatrix;
      GenerateMassMatrix::computeMassMatrix(precomputationState.simulationMesh, &massMatrix, true); 

      if (!fullMatrix)
      {
        // constrain the degrees of freedom
        int numConstrainedVertices = (int) (precomputationState.fixedVertices.size());
        int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
        int i = 0;
        for(set<int> :: iterator iter = precomputationState.fixedVertices.begin(); 
          iter != precomputationState.fixedVertices.end(); iter++)
        {
          constrainedDOFs[3*i+0] = 3 * (*iter) + 1;
          constrainedDOFs[3*i+1] = 3 * (*iter) + 2;
          constrainedDOFs[3*i+2] = 3 * (*iter) + 3;
          i++;
        }

        int oneIndexed = 1;
        massMatrix->RemoveRowsColumns(
          3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
        free(constrainedDOFs);
      }

      const char * filename = massMatrixFilename.mb_str();
      int code = massMatrix->Save((char*)filename);
      delete(massMatrix);
      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save mass matrix to ") + massMatrixFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnDisplayRigidModes(wxCommandEvent& event)
{
  string information = string("Rigid modes: ") + string("\n"); 

  if (precomputationState.numRigidModes == 0)
    information = information + "none\n";
  else
    for(int i=0; i<precomputationState.numRigidModes; i++)
      information = information + int2string(i+1) + string("\n");

  wxWindow * parent = NULL;
  WxDialogWrapper<wxMessageDialog> dlg( new wxMessageDialog(parent,
                wxString(information.c_str(), wxConvUTF8),
                wxString(_T("Rigid mode information"), wxConvUTF8),
                wxOK | wxICON_INFORMATION) );
  dlg.get()->ShowModal();
}

void MyFrame::OnSetRigidModes(wxCommandEvent& event)
{
  wxDialog * dlg = new wxDialog(this, -1, _T("Select number of rigid modes"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );

  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * numModesText = new wxStaticText(dlg, -1,
       _T("Number of rigid modes: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  char numRigidModesStringC[96];
  sprintf(numRigidModesStringC, "%d", precomputationState.numRigidModes);
  wxTextCtrl * numModesControl = new wxTextCtrl(dlg, -1,
      wxString(numRigidModesStringC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numModesSizer = new wxBoxSizer(wxHORIZONTAL);
  numModesSizer->Add(numModesText, 0, wxCENTER);
  numModesSizer->Add(numModesControl, 0, wxCENTER);
  dlgSizer->Add(numModesSizer, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString valueString = numModesControl->GetValue();

  delete(dlg);

  long value;
  bool goodInput = valueString.ToLong(&value);

  if (goodInput)
  {
    if ((value < 0) || (value > 6))
      goodInput = false;
  }

  if (goodInput)
  {
    precomputationState.numRigidModes = value;
  }
  else
  {
    this->errMsg( _T("Invalid number of rigid modes"),  _T("Invalid number of rigid modes: ") +  valueString );
  }
}

