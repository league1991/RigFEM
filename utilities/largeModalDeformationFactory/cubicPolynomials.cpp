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

// cubic polynomial computation

#include "StVKReducedInternalForces.h"
#include "StVKReducedInternalForcesWX.h"
#include "StVKCubeABCD.h"
#include "matrixIO.h"
#include "StVKElementABCDLoader.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnLoadCubicPolynomials(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load cubic polynomials"),
    uiState.currentWorkingDirectory, _T(""), _T("Cubic polynomials files(*.cub)|*.cub|All files(*.*)|*.*"),
    wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString cubicPolynomialsFilename = dlg->GetPath();
    SaveCurrentWorkingDirectory(cubicPolynomialsFilename);
    if( !cubicPolynomialsFilename.empty() )
    {
      StVKReducedInternalForces * newCubicPolynomials = NULL;

      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = cubicPolynomialsFilename.mb_str();
      newCubicPolynomials = new StVKReducedInternalForces((char*) filename);
      int code = 0; // should also catch any errors
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  
        _T("Unable to load cubic polynomials from ") + cubicPolynomialsFilename );
        dlg->Destroy();
        return;
      }

      // success
      delete(precomputationState.cubicPolynomials);
      precomputationState.cubicPolynomials = newCubicPolynomials;

      precomputationState.cubicPolynomialsAvailable = true;

      UpdateMenus();

      Refresh();
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSaveCubicPolynomials(wxCommandEvent& event)
{
   wxFileDialog *dlg = new wxFileDialog(this, _T("Save cubic polynomials"),
     uiState.currentWorkingDirectory, _T(""), _T("Cubic Polynomials Files(*.cub)|*.cub|All files(*.*)|*.*"),
	wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);
   if ( dlg->ShowModal() == wxID_OK )
   {
     wxString cubicPolynomialsFilename( dlg->GetPath() );
     SaveCurrentWorkingDirectory(cubicPolynomialsFilename);
     if( !cubicPolynomialsFilename.empty() )
     {
       const char * filename = cubicPolynomialsFilename.mb_str();
       int code = precomputationState.cubicPolynomials->Save((char*)filename);

       if (code != 0)
       {
         this->errMsg( _T("Saving error"),  
           _T("Unable to save cubic polynomials to ") + cubicPolynomialsFilename );
         dlg->Destroy();
         return;
       }	
     }
   }

   dlg->Destroy();
}

void MyFrame::OnComputeCubicPolynomials(wxCommandEvent& event)
{
  if (!precomputationState.nonLinearModesAvailable)
  {
    this->errMsg( _T("Cannot compute cubic polynomials"),  
      _T("Nonlinear modes are not available.\n"
      ) );
    return;
  }

  wxDialog * dlg = new wxDialog(this, -1, _T("Select the number of threads"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );

  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);
  wxStaticText * numComputationThreadsText = new wxStaticText(dlg, -1,
       _T("Number of computation threads: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  char numComputationThreadsC[96];
  sprintf(numComputationThreadsC, "%d", uiState.numComputationThreads);
  wxTextCtrl * numComputationThreadsControl = new wxTextCtrl(dlg, -1,
      wxString(numComputationThreadsC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numComputationThreadsSizer = new wxBoxSizer(wxHORIZONTAL);
  numComputationThreadsSizer->Add(numComputationThreadsText, 0, wxCENTER);
  numComputationThreadsSizer->Add(numComputationThreadsControl, 0, wxCENTER);
  dlgSizer->Add(numComputationThreadsSizer, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString valueString = numComputationThreadsControl->GetValue();

  delete(dlg);

  long value;
  bool goodInput = valueString.ToLong(&value);

  if (goodInput)
  {
    if ((value <= 0) || (value > 65536))
      goodInput = false;
  }

  if (goodInput)
  {
    uiState.numComputationThreads = value;

    StVKReducedInternalForces * newCubicPolynomials = NULL;

    SetCursor(*wxHOURGLASS_CURSOR);
    int code;
    CubicPolynomialsWorker(&code, &newCubicPolynomials);
    SetCursor(*wxSTANDARD_CURSOR);

    if (code != 0)
    {
      this->errMsg( _T("Cubic polynomials computation failed"),  
             _T("Cubic polynomials computation failed.") );
      free(newCubicPolynomials);
      return;
    }
    else
    {
      delete(precomputationState.cubicPolynomials);
      precomputationState.cubicPolynomials = newCubicPolynomials;
      precomputationState.cubicPolynomialsAvailable = true;

      UpdateMenus();

      Refresh();
    }
  }
}

void * MyFrame::CubicPolynomialsWorker(int * code, StVKReducedInternalForces ** newCubicPolynomial)
{
  // compute cubic polynomials
  StVKElementABCD * precomputedABCDIntegrals = StVKElementABCDLoader::load(precomputationState.simulationMesh);

  if (uiState.numComputationThreads == 1)
  {
    *newCubicPolynomial = new StVKReducedInternalForces(
      precomputationState.rNonLin,
      precomputationState.nonLinearModalMatrix->GetMatrix(),
      precomputationState.simulationMesh, precomputedABCDIntegrals); 
  }
  else
  {
    *newCubicPolynomial = new StVKReducedInternalForcesWX(
      precomputationState.rNonLin,
      precomputationState.nonLinearModalMatrix->GetMatrix(),
      precomputationState.simulationMesh, precomputedABCDIntegrals, false, uiState.numComputationThreads); 
  }

  delete(precomputedABCDIntegrals);

  computationRunning = -1;

  *code = 0;
  return NULL;
}

