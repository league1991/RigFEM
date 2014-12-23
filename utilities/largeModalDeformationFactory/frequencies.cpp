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

// load/save/scale frequencies

#include <iomanip>
#include <sstream>
#include <iostream>
#include "sparseMatrix.h"
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "matrixIO.h"
#include "insertRows.h"
#include "volumetricMeshENuMaterial.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnLoadFrequencies(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load frequencies"), uiState.currentWorkingDirectory, _T(""), _T("Frequency Files(*.freq)|*.freq|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString frequencyFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(frequencyFilename);
    if( !frequencyFilename.empty() )
    {
      int newr;
      double * newFrequencies = NULL;

      int m1;
      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = frequencyFilename.mb_str();
      int code = ReadMatrixFromDisk((char*)filename, &newr, &m1, &newFrequencies);
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  _T("Unable to load frequencies from ") + frequencyFilename );
        dlg->Destroy();
        return;
      }

      if (precomputationState.linearModesAvailable && (newr != precomputationState.rLin))
      {
        char s[4096];
        const char * filename = frequencyFilename.mb_str();
        sprintf(s, "The number of frequencies (%d) "
          "in %s does not match the currently available number of linear modes (%d).",
          newr, (char*)filename, precomputationState.rLin);
        this->errMsg( _T("Loading error"), wxString(s, wxConvUTF8));
        free(newFrequencies);
        dlg->Destroy();
        return;
      }

      // success
      free(precomputationState.frequencies);
      precomputationState.rLin = newr;
      precomputationState.frequencies = newFrequencies;

      precomputationState.frequenciesAvailable = true;
      UpdateMenus();
      UpdateModalToolbar();

      Refresh();
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSaveFrequencies(wxCommandEvent& event)
{
   wxFileDialog *dlg = new wxFileDialog(this, _T("Save frequencies"), uiState.currentWorkingDirectory, _T(""), _T("Frequency Files(*.freq)|*.freq|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString frequencyFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(frequencyFilename);
    if( !frequencyFilename.empty() )
    {
      const char * filename = frequencyFilename.mb_str();
      int code = WriteMatrixToDisk((char*)filename,
        precomputationState.rLin, 1, precomputationState.frequencies);

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save frequencies to ") + frequencyFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnExportFrequencies(wxCommandEvent& event)
{
   wxFileDialog *dlg = new wxFileDialog(this, _T("Export frequencies"), uiState.currentWorkingDirectory, _T(""), _T("Text Files(*.txt)|*.txt|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString frequencyFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(frequencyFilename);
    if( !frequencyFilename.empty() )
    {
      const char * filename = frequencyFilename.mb_str();
      FILE * fout = fopen((char*)filename, "w");
      if (fout)
      {
        fprintf(fout,"%d\n", precomputationState.rLin);
        for(int j=0; j<precomputationState.rLin; j++)
        {
          fprintf(fout, "%.15f\n", precomputationState.frequencies[j]);
        }
        fclose(fout);
      }
      else
      {
        this->errMsg( _T("Exporting error"),  
          _T("Unable to export frequencies to ") + frequencyFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnScaleTo1Hz(wxCommandEvent& event)
{
  wxMessageDialog * confirmationDialog = new wxMessageDialog
    (this, _T("Scale frequency spectrum so that lowest frequency is at 1.0 Hz ?"),
      _T("Scale frequency spectrum"), 
      wxYES_NO | wxICON_QUESTION);

  if (confirmationDialog->ShowModal() == wxID_YES)
  {
    double factor = 1.0 / precomputationState.frequencies[precomputationState.numRigidModes];
    for(int i=0; i<precomputationState.rLin; i++)
      precomputationState.frequencies[i] *= factor;

    ScaleYoungsModulus(factor * factor);
  }

  UpdateModalToolbar();

  delete(confirmationDialog);
}

void MyFrame::OnScaleFrequencies(wxCommandEvent& event)
{
  char oldFactorStringC[4096];
  sprintf(oldFactorStringC,"%f", uiState.lastFrequencyScaleFactor);

  wxDialog * dlg = new wxDialog(this, -1, _T("Scale frequency spectrum"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * factorText = new wxStaticText(dlg, -1, 
       _T("Multiplicative factor:"), wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * factorControl = new wxTextCtrl(dlg, -1, 
      wxString(oldFactorStringC, wxConvUTF8), wxDefaultPosition, wxSize(80,-1), 
      0,  wxDefaultValidator);

  wxBoxSizer * factorControlSizer = new wxBoxSizer(wxHORIZONTAL);
  factorControlSizer->Add(factorText, 0, wxCENTER);
  factorControlSizer->Add(factorControl, 0, wxCENTER);
  dlgSizer->Add(factorControlSizer, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, 
    wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  double multiplicativeFactor;

  if (!(factorControl->GetValue().ToDouble(&multiplicativeFactor)))
  {
    errMsg(_T("Invalid scaling factor"), _T("Invalid scaling factor: ") + factorControl->GetValue());
    delete(dlg);
    return;
  }

  delete(dlg);

  for(int i=0; i<precomputationState.rLin; i++)
    precomputationState.frequencies[i] *= multiplicativeFactor;

  ScaleYoungsModulus(multiplicativeFactor * multiplicativeFactor);

  uiState.lastFrequencyScaleFactor = multiplicativeFactor;

  UpdateModalToolbar();
}

void MyFrame::ScaleYoungsModulus(double factor)
{
  for(int i=0; i < precomputationState.simulationMesh->getNumMaterials(); i++)
  {
    VolumetricMesh::Material * material = precomputationState.simulationMesh->getMaterial(i);
    VolumetricMesh::ENuMaterial * eNuMaterial =  downcastENuMaterial(material);
    if (eNuMaterial == NULL)
      printf("Internal error: eNuMaterial is NULL);");
    else
      eNuMaterial->setE(factor * eNuMaterial->getE());
  }
}

string MyFrame::double2string(double d)
{
  ostringstream output;
  output << setprecision(6) << d;
  return output.str();
}

string MyFrame::int2string(int n)
{
  ostringstream output;
  output << n;
  return output.str();
}

void MyFrame::OnDisplayFrequencies(wxCommandEvent& event)
{
  string information = string("Num frequencies: ") + int2string(precomputationState.rLin) + string("\n\n");

  if (precomputationState.rLin <= 20)
  {
    for(int i=0; i<precomputationState.rLin; i++)
      information = information + double2string(precomputationState.frequencies[i]) + string("\n");
  }
  else
  {
    for(int i=0; i<18; i++)
    {
      information = information + double2string(precomputationState.frequencies[i]) + string("\n");
    }
    information = information + string("...\n");
    information = information + double2string(precomputationState.frequencies[precomputationState.rLin-2]) + string("\n");
    information = information + double2string(precomputationState.frequencies[precomputationState.rLin-1]) + string("\n");
  }

  wxWindow * parent = NULL;
  WxDialogWrapper<wxMessageDialog> dlg( new wxMessageDialog(parent,
    wxString(information.c_str(), wxConvUTF8), 
    wxString(_T("Frequency information"), wxConvUTF8),
    wxOK | wxICON_INFORMATION) );
  dlg.get()->ShowModal();
}

void MyFrame::OnDeleteFrequencies(wxCommandEvent& event)
{
  char eraseRangeLoC[96];
  char eraseRangeHiC[96];
  sprintf(eraseRangeLoC, "%d", uiState.eraseRangeLo+1);
  sprintf(eraseRangeHiC, "%d", uiState.eraseRangeHi);

  wxDialog * dlg = new wxDialog(this, -1, _T("Select the frequency range to be erased"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );

  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * staticText1 = new wxStaticText(dlg, -1,
       _T("Erase frequencies from: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxStaticText * staticText2 = new wxStaticText(dlg, -1,
       _T(" to: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * eraseRangeLoControl = new wxTextCtrl(dlg, -1,
      wxString(eraseRangeLoC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxTextCtrl * eraseRangeHiControl = new wxTextCtrl(dlg, -1,
      wxString(eraseRangeHiC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * sizer1 = new wxBoxSizer(wxHORIZONTAL);
  sizer1->Add(staticText1, 0, wxCENTER);
  sizer1->Add(eraseRangeLoControl, 0, wxCENTER);
  sizer1->Add(staticText2, 0, wxCENTER);
  sizer1->Add(eraseRangeHiControl, 0, wxCENTER);
  dlgSizer->Add(sizer1, 0, wxALL | wxCENTER, 8);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString valueStringLo = eraseRangeLoControl->GetValue();
  wxString valueStringHi = eraseRangeHiControl->GetValue();

  delete(dlg);

  long valueLo;
  bool goodInput = valueStringLo.ToLong(&valueLo);
  long valueHi;
  goodInput = goodInput && valueStringHi.ToLong(&valueHi);

  if (goodInput)
  {
    if ((valueLo <= 0) || (valueLo > precomputationState.rLin))
      goodInput = false;

    if ((valueHi <= 0) || (valueHi > precomputationState.rLin))
      goodInput = false;

    if (valueLo > valueHi)
      goodInput = false;
  }

  if (!goodInput)
  {
    this->errMsg( _T("Incorrect frequency range"),
          _T("Error: Incorrect frequency range"));
    return;
  }

  uiState.eraseRangeLo = valueLo - 1;
  uiState.eraseRangeHi = valueHi;

  int numErasedFrequencies = uiState.eraseRangeHi - uiState.eraseRangeLo;
  int n3 = 3 * precomputationState.simulationMesh->getNumVertices();

  double * U = (double*) malloc (sizeof(double) * n3 * precomputationState.rLin);
  memcpy(U, (precomputationState.linearModalMatrix)->GetMatrix(), sizeof(double) * n3 * precomputationState.rLin);

  // erase frequencies and linear modes
  for(int i=uiState.eraseRangeHi; i<precomputationState.rLin; i++)
  {
    precomputationState.frequencies[i-numErasedFrequencies] = precomputationState.frequencies[i];
    memcpy(&U[ELT(n3,0,i-numErasedFrequencies)], &U[ELT(n3,0,i)], sizeof(double) * n3);
  }

  precomputationState.rLin -= numErasedFrequencies;

  precomputationState.frequencies = (double*) realloc (precomputationState.frequencies, sizeof(double) * precomputationState.rLin);

  delete(precomputationState.linearModalMatrix);
  int noInternalCopy = 1;
  precomputationState.linearModalMatrix = new ModalMatrix(n3 / 3, precomputationState.rLin, U, noInternalCopy);

  if (modeSelectionControl->GetValue() > precomputationState.rLin)
    modeSelectionControl->SetValue(precomputationState.rLin);
  modeSelectionControl->SetRange(1, precomputationState.rLin);

  UpdateModalToolbar();
}

