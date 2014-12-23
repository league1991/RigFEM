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

// interpolation from volumetric to rendering mesh

#include "matrixIO.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnInterpolateLinearModes(wxCommandEvent& event)
{
  InterpolateMatrix(_T("linear modes"), precomputationState.linearModalMatrix);
}

void MyFrame::OnInterpolateNonLinearModes(wxCommandEvent& event)
{
  InterpolateMatrix(_T("nonlinear modes"), precomputationState.nonLinearModalMatrix);
}

void MyFrame::OnInterpolateModalDerivatives(wxCommandEvent& event)
{
  InterpolateMatrix(_T("modal derivatives"), precomputationState.modalDerivativesMatrix);
}

void MyFrame::OnInterpolateSimulationData(wxCommandEvent& event)
{
  InterpolateMatrix(_T("external simulation data"), precomputationState.sketchDataMatrix);
}

void MyFrame::InterpolateMatrix(wxString dataDescription, ModalMatrix * inputModalMatrix)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Interpolate ") + dataDescription + _T(" to triangle mesh"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.U)|*.U|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString outputFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(outputFilename);
    if( !outputFilename.empty() )
    {
      SetCursor(*wxHOURGLASS_CURSOR);
      if (!precomputationState.interpolationDataAvailable)
        BuildInterpolant();

      // interpolate
      printf("Interpolating data...\n");
      double * inputMatrix = inputModalMatrix->GetMatrix();
      int nTarget = (int)(precomputationState.renderingMesh->getNumVertices());
      double * outputMatrix = (double*) malloc (sizeof(double) * 3 * nTarget * inputModalMatrix->Getr());
      for(int i=0; i<inputModalMatrix->Getr(); i++)
      {
        precomputationState.simulationMesh->interpolate(
          &inputMatrix[ELT(3*precomputationState.simulationMesh->getNumVertices(),0,i)],
          &outputMatrix[ELT(3*nTarget,0,i)],
          nTarget, precomputationState.simulationMesh->getNumElementVertices(),
          precomputationState.interpolationData_vertices,
          precomputationState.interpolationData_weights);
      }

      // save file to disk
      const char * filename = outputFilename.mb_str();
	  printf("Saving output to %s.\n", (char*)filename);
	  int code = WriteMatrixToDisk((char*)filename,
        3 * nTarget, 
        inputModalMatrix->Getr(), 
        outputMatrix);

      free(outputMatrix);

      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save interpolated data to ") + outputFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::BuildInterpolant()
{
  printf("Building interpolation data structure...\n");
  int nTarget;
  double * targetLocations;
  precomputationState.renderingMesh->exportGeometry
    (&nTarget, &targetLocations, NULL, NULL);

  int numExternalVertices = 
    precomputationState.simulationMesh->generateInterpolationWeights(
        nTarget, targetLocations, 
        &(precomputationState.interpolationData_vertices),
        &(precomputationState.interpolationData_weights),
	-1, NULL, 1);

  printf("Detected %d rendering mesh vertices outside of the simulation mesh.\n", 
    numExternalVertices);
          
  free(targetLocations);
}

void MyFrame::OnSaveInterpolant(wxCommandEvent& event)
{
  wxFileDialog * dlg = new wxFileDialog(this, _T("Save interpolant to file"), uiState.currentWorkingDirectory, _T(""), _T("Interpolant Text Files(*.interp)|*.interp|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString outputFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(outputFilename);
    if( !outputFilename.empty() )
    {
      SetCursor(*wxHOURGLASS_CURSOR);
      if (!precomputationState.interpolationDataAvailable)
        BuildInterpolant();

      int nTarget = (int)(precomputationState.renderingMesh->getNumVertices());
      const char * filename = outputFilename.mb_str();

      int code = precomputationState.simulationMesh->saveInterpolationWeights
        ((char*)filename, nTarget, 
         precomputationState.simulationMesh->getNumElementVertices(),
         precomputationState.interpolationData_vertices, 
         precomputationState.interpolationData_weights);

      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save interpolated data to ") + outputFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();  
}

