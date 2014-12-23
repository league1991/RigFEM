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

// nonlinear mode computation

#include "sparseMatrix.h"
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "StVKHessianTensor.h"
#include "insertRows.h"
#include "sparseSolvers.h"
#include "matrixIO.h"
#include "matrixPCA.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnLoadNonLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load nonlinear modes"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.U)|*.U|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString nonLinearModesFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(nonLinearModesFilename);
    if( !nonLinearModesFilename.empty() )
    {
      int newrNonLin;
      double * newNonLinearModes = NULL;

      int n1;
      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = nonLinearModesFilename.mb_str();
      int code = ReadMatrixFromDisk((char*)filename, &n1, &newrNonLin, &newNonLinearModes);
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  
          _T("Unable to load nonlinear modes from ") + nonLinearModesFilename );
        dlg->Destroy();
        return;
      }

      if (n1 != 3 * precomputationState.simulationMesh->getNumVertices())
      {
        this->errMsg( _T("Loading error"),  
          _T("The number of vertices in ") + nonLinearModesFilename + _T(" does not match the simulation mesh."));
        free(newNonLinearModes);
        dlg->Destroy();
        return;
      }

      // success
      delete(precomputationState.nonLinearModalMatrix);
      precomputationState.rNonLin = newrNonLin;
      precomputationState.nonLinearModalMatrix = new ModalMatrix(precomputationState.simulationMesh->getNumVertices(), precomputationState.rNonLin, newNonLinearModes);
      free(newNonLinearModes);

      precomputationState.nonLinearModesAvailable = true;

      modeSelectionControl->SetValue(1);

      SelectView(UIState::VIEW_NONLINEAR_MODES);
      SetAutoRenderingMagnitude(precomputationState.nonLinearModalMatrix);

      UpdateMenus();

      myGLCanvas->UpdateNonLinearModesRenderData();

      Refresh();
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSaveNonLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Save nonlinear modes"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.U)|*.U|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString nonLinearModesFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(nonLinearModesFilename);
    if( !nonLinearModesFilename.empty() )
    {
      const char * filename = nonLinearModesFilename.mb_str();
      int code = WriteMatrixToDisk((char*)filename,
        3 * precomputationState.nonLinearModalMatrix->Getn(), 
        precomputationState.nonLinearModalMatrix->Getr(), 
        precomputationState.nonLinearModalMatrix->GetMatrix());

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
            _T("Unable to save nonlinear modes to ") + nonLinearModesFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnExportNonLinearModes(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Export linear modes"), uiState.currentWorkingDirectory, _T(""), _T("Text Files(*.txt)|*.txt|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString nonLinearModesFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(nonLinearModesFilename);
    if( !nonLinearModesFilename.empty() )
    {
      const char * filename = nonLinearModesFilename.mb_str();
      FILE * fout = fopen((char*)filename, "w");
      if (fout)
      {
        double * U = precomputationState.nonLinearModalMatrix->GetMatrix();
        int n = precomputationState.nonLinearModalMatrix->Getn();
        fprintf(fout,"%d\n%d\n", 3*n, precomputationState.rNonLin);
        for(int i=0; i<3*n; i++)
        {
          for(int j=0; j<precomputationState.rNonLin; j++)
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
         _T("Unable to export nonlinear modes to ") + nonLinearModesFilename );
       dlg->Destroy();
       return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnComputeNonLinearModes(wxCommandEvent& event)
{
  if (! ( (precomputationState.linearModesAvailable && 
           precomputationState.frequenciesAvailable &&
           precomputationState.modalDerivativesAvailable)
         || (precomputationState.sketchDataAvailable)
         )
     )
  {
    this->errMsg( _T("Cannot compute the simulation basis"),  
      _T("You must compute/provide at least one of the following:\n"
         L"  (1) linear modes, frequencies, and modal derivatives, OR\n"
         L"  (2) external simulation data.\n"
      ) );
      return;
  }

  char numNonLinearModesStringC[256];
  sprintf(numNonLinearModesStringC, "%d", uiState.numComputedNonLinearModes);

  wxDialog * dlg = new wxDialog(this, -1, _T("Compute nonlinear modes (via SVD)"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  // create dialog box to choose source of data
  wxString radioBoxChoices[2] = { 
    _T("Linear modes and modal derivatives"), _T("External simulation data") };

  wxRadioBox * dataOriginRadioBox = new wxRadioBox(dlg, -1, 
    _T("Source of data for nonlinear modes:"),
    wxDefaultPosition, wxDefaultSize, 2, radioBoxChoices, 
    2, wxRA_SPECIFY_ROWS);

  dataOriginRadioBox->Enable(0, precomputationState.modalDerivativesAvailable);
  dataOriginRadioBox->Enable(1, precomputationState.sketchDataAvailable);

  if (precomputationState.modalDerivativesAvailable)
    dataOriginRadioBox->SetSelection(0);
  else
    dataOriginRadioBox->SetSelection(1);

  dlgSizer->Add(dataOriginRadioBox, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  wxStaticText * numModesText = new wxStaticText(dlg, -1, 
       _T("Number of nonlinear modes: "), wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numModesControl = new wxTextCtrl(dlg, -1, 
      wxString(numNonLinearModesStringC, wxConvUTF8), 
      wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numModesSizer = new wxBoxSizer(wxHORIZONTAL);
  numModesSizer->Add(numModesText, 0, wxCENTER);
  numModesSizer->Add(numModesControl, 0, wxCENTER);
  dlgSizer->Add(numModesSizer, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, 
    wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString valueString = numModesControl->GetValue();
  int dataOrigin = dataOriginRadioBox->GetSelection();

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
    uiState.numComputedNonLinearModes = value;

    double * newNonLinearModes = NULL;

    SetCursor(*wxHOURGLASS_CURSOR);
    int code;
    NonLinearModesWorker(&code, dataOrigin,
      uiState.numComputedNonLinearModes,
      &newNonLinearModes );
    SetCursor(*wxSTANDARD_CURSOR);

    if (code != 0)
    {
      char s[96];
      if (code < 0)
        sprintf(s, "Nonlinear mode computation failed. Memory allocation error.");
      else
        sprintf(s, "Nonlinear mode computation failed. dgesvd exit code: %d\n", code);
      this->errMsg( _T("Nonlinear mode computation failed"),  wxString(s, wxConvUTF8));
      free(newNonLinearModes);
      return;
    }
    else
    {
      precomputationState.nonLinearModesAvailable = true;
      UpdateMenus();
 
      delete(precomputationState.nonLinearModalMatrix);
      precomputationState.rNonLin = uiState.numComputedNonLinearModes;
      precomputationState.nonLinearModalMatrix = new ModalMatrix(
        precomputationState.simulationMesh->getNumVertices(), 
        precomputationState.rNonLin, newNonLinearModes);
      free(newNonLinearModes);

      modeSelectionControl->SetValue(1);

      SelectView(UIState::VIEW_NONLINEAR_MODES);
      SetAutoRenderingMagnitude(precomputationState.nonLinearModalMatrix);

      UpdateMenus();

      myGLCanvas->UpdateNonLinearModesRenderData();

      Refresh();
    }
  }
  else
  {
    this->errMsg( _T("Invalid number of nonlinear modes"),  _T("Invalid number of nonlinear modes: ") + valueString );
  }
}

void * MyFrame::NonLinearModesWorker(int * code, int dataOrigin, int numNonLinearModes, double ** modes_ )
{
  int n3 = 3 * precomputationState.simulationMesh->getNumVertices();
  int numDataVectors = 0;

  // compute the lumped mass matrix
  SparseMatrix * massMatrix; // will be sparse n3 x n3
  GenerateMassMatrix::computeMassMatrix( precomputationState.simulationMesh, &massMatrix, true); // exitCode will always be 0

  // mass-normalize modal derivatives
  double * modalDerivatives = precomputationState.modalDerivativesMatrix->GetMatrix();
  double * normalizedModalDerivatives = (double*) malloc (sizeof(double) * n3 * precomputationState.numDeriv);
  memcpy(normalizedModalDerivatives, modalDerivatives, sizeof(double) * n3 * precomputationState.numDeriv);

  for(int i=0; i < precomputationState.numDeriv; i++)
    massMatrix->NormalizeVector(&(normalizedModalDerivatives[n3 * i]));

  double * dataMatrix;
  if (dataOrigin == 0)
  {
    // use linear modes and derivatives
    // construct PCA data matrix:
    int numUsedLinearModes = precomputationState.rLin - precomputationState.numRigidModes;
    numDataVectors = numUsedLinearModes + numUsedLinearModes * (numUsedLinearModes + 1) / 2;
    printf("Number of PCA datamatrix columns: %d.\n", numDataVectors);
    dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);

    printf("Generating datamatrix for SVD...\n");

    double lambda0 = precomputationState.frequencies[precomputationState.numRigidModes] * precomputationState.frequencies[precomputationState.numRigidModes];

    // scale linear modes 
    double * Ulin = precomputationState.linearModalMatrix->GetMatrix();
    for(int i=0; i<numUsedLinearModes; i++)
    {
      double lambda = precomputationState.frequencies[precomputationState.numRigidModes + i] * precomputationState.frequencies[precomputationState.numRigidModes + i];
      double factor = lambda0 / lambda;
      for(int vertex=0; vertex< n3; vertex++)
        dataMatrix[ELT(n3,vertex,i)] = factor * Ulin[ELT(n3,vertex,precomputationState.numRigidModes + i)];
    } 

    // scale modal derivatives
    int pos = 0;
    for(int i=0; i<numUsedLinearModes; i++)
    {
      for(int j=i; j<numUsedLinearModes; j++)
      {
        double lambdai = precomputationState.frequencies[precomputationState.numRigidModes + i] * precomputationState.frequencies[precomputationState.numRigidModes + i];
        double lambdaj = precomputationState.frequencies[precomputationState.numRigidModes + j] * precomputationState.frequencies[precomputationState.numRigidModes + j];
        double factor = lambda0 * lambda0 / (lambdai * lambdaj);

        for(int vertex=0; vertex < n3; vertex++)
          dataMatrix[ELT(n3, vertex, numUsedLinearModes + pos)] = 
            factor * normalizedModalDerivatives[ELT(n3, vertex, pos)];
        pos++;
      }
    }
  }
  else
  {
    // data from external simulation
    numDataVectors = precomputationState.sketchDataMatrix->Getr();
    dataMatrix = (double*) malloc (sizeof(double) * n3 * numDataVectors);
    memcpy(dataMatrix, precomputationState.sketchDataMatrix->GetMatrix(), sizeof(double) * n3 * numDataVectors);
  }

  free(normalizedModalDerivatives);

  // do lumped-mass-PCA on dataMatrix ( n3 x numDataVectors )
  
  double * ones = (double*) malloc (sizeof(double) * n3);
  for(int i=0; i<n3; i++)
    ones[i] = 1.0;

  double * LTDiagonal = (double*) malloc (sizeof(double) * n3);
  massMatrix->MultiplyVector(ones, LTDiagonal);
  free(ones);
  delete(massMatrix);

  // sqrt
  for(int i=0; i<n3; i++)
    LTDiagonal[i] = sqrt(LTDiagonal[i]);

  // number of retained dimensions can't be more than num linear modes + num derivatives
  if (uiState.numComputedNonLinearModes > numDataVectors)
    uiState.numComputedNonLinearModes = numDataVectors;

  // premultiply by LT
  for(int i=0; i<n3; i++)
    for(int j=0; j < numDataVectors; j++)
      dataMatrix[ELT(n3, i, j)] *= LTDiagonal[i];

  // do SVD on dataMatrix ( n3 x numDataVectors ), retain uiState.numComputedNonLinearModes modes

  ThresholdingSpecification thresholdingSpecification;
  thresholdingSpecification.tresholdingType = ThresholdingSpecification::numberOfModesBased;
  thresholdingSpecification.rDesired = uiState.numComputedNonLinearModes;

  int outputr;
  int matrixPCACode = 0;
  if ( ((matrixPCACode = MatrixPCA(
    &thresholdingSpecification, n3, numDataVectors, dataMatrix, &outputr)) != 0) 
    || (outputr != uiState.numComputedNonLinearModes))
  {
    printf("Error performing SVD. Code: %d\n", matrixPCACode);
    *code = matrixPCACode;
    free(dataMatrix);
    free(LTDiagonal);
    return NULL;
  }

  // solve L^T U = V
  for(int i=0; i<n3; i++)
    for(int j=0; j < uiState.numComputedNonLinearModes; j++)
      dataMatrix[ELT(n3, i, j)] /= LTDiagonal[i];

  free(LTDiagonal);

  // export data
  *modes_ = (double*) realloc (dataMatrix, 
    sizeof(double) * n3 * uiState.numComputedNonLinearModes);

  computationRunning = -1;

  *code = 0;
  return NULL;
}

