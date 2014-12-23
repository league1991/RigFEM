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

// create the configuration files for a runtime simulation

#include "StVKReducedInternalForces.h"
#include "StVKReducedStiffnessMatrix.h"
#include "matrixIO.h"
#include "matrix.h"
#include "generateMassMatrix.h"
#include "largeModalDeformationFactory.h"

static char * folder_open_xpm_2[] = {
/* columns rows colors chars-per-pixel */
"16 15 31 1",
"6 c #9BACC2",
"w c #547B99",
"5 c #94A5BD",
". c #376485",
"; c #F1F4F7",
"o c #7F99B4",
"2 c #D1D9E5",
"- c #EAEDF3",
"O c #718BA7",
"0 c #65839D",
"* c #DCE2EA",
": c #F5F6F7",
"7 c #597B9A",
"X c #8DA0B9",
"  c None",
"+ c #467291",
"q c #305F81",
"& c #D6DFE7",
"3 c #6A89A2",
"1 c #A8B6CA",
"= c #E4E9ED",
"> c #F8F9FA",
", c #FDFDFE",
"9 c #215579",
"8 c #7F97B0",
"@ c #B3BFD1",
"< c #7A90AC",
"$ c #C2CBDB",
"4 c #A2B3C5",
"% c #CAD6E1",
"# c #BBC4D6",
/* pixels */
"                L",
".....           L",
".XXXo.          L",
".XXXXO........  L",
".XXXXXXXXXXXX.  L",
".XXXXXXXXXXXX.  L",
".X++++++++++++++",
".X+@#$%&*=-;:>,+",
".<.1@#$%2*=-;:23",
"..X41@#$%2*=-;3 ",
"..X561@#$%2*=-3 ",
".78X561@#$%2*%3 ",
"90<8X561@#$%23  L",
"q++++++++++++w  L",
"                L"
};


const char lightingString[] = 
"*globalAmbientIntensity\n"
"0.2\n"
"\n"
"*enableSpecularTerm\n"
"false\n"
"\n"
"#### light 3\n"
"\n"
"*lightEnabled_3\n"
"true\n"
"\n"
"*position_3_X\n"
"-30\n"
"*position_3_Y\n"
"30\n"
"*position_3_Z\n"
"-30\n"
"\n"
"*lightIntensity_3\n"
"0.5\n"
"\n"
"#### light 6\n"
"\n"
"*lightEnabled_6\n"
"true\n"
"\n"
"*position_6_X\n"
"30\n"
"*position_6_Y\n"
"30\n"
"*position_6_Z\n"
"30\n"
"\n"
"*lightIntensity_6\n"
"0.5\n";

int MyFrame::PreprocessSingleSimulation(wxString & projectName, int numLinearModes, int numNonLinearModes, int numComputationThreads)
{
  // === compute linear modes

  int oldNumLinearModes = uiState.numComputedLinearModes;
  uiState.numComputedLinearModes = numLinearModes;

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
    uiState.numComputedLinearModes = oldNumLinearModes;
    this->errMsg( _T("Linear mode computation failed"),
           _T("Linear mode computation failed.") );
    free(newFrequencies);
    free(newLinearModes);
    return 1;
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

    //uiState.numComputedNonLinearModes = 2 * (precomputationState.rLin - precomputationState.numRigidModes);

    modeSelectionControl->SetValue(1);

    SelectView(UIState::VIEW_LINEAR_MODES);
    SetAutoRenderingMagnitude(precomputationState.linearModalMatrix);

    UpdateMenus();

    myGLCanvas->UpdateLinearModesRenderData();

    Refresh();
  }

  // === compute modal derivatives

  double * modalDerivatives = NULL;
  int code;
  SetCursor(*wxHOURGLASS_CURSOR);
  ComputeModalDerivatives(&code, &modalDerivatives);
  SetCursor(*wxSTANDARD_CURSOR);

  if (code != 0)
  {
    this->errMsg( _T("Modal derivatives computation failed"),
           _T("Modal derivatives computation failed.") );
    free(modalDerivatives);
    return 2;
  }
  else
  {
    printf("Modal derivative computation succeeded.\n"); fflush(NULL);

    precomputationState.modalDerivativesAvailable = true;

    delete(precomputationState.modalDerivativesMatrix);
    precomputationState.modalDerivativesMatrix = new ModalMatrix(
        precomputationState.simulationMesh->getNumVertices(),
        (precomputationState.rLin - precomputationState.numRigidModes) * (precomputationState.rLin - precomputationState.numRigidModes + 1) / 2,
        modalDerivatives);
    free(modalDerivatives);

    UpdateMenus();
  }

  // === compute nonlinear modes 

  int oldNumNonLinearModes = uiState.numComputedNonLinearModes;
  uiState.numComputedNonLinearModes = numNonLinearModes;

  double * newNonLinearModes = NULL;

  int dataOrigin = 0; // linear modes and modal derivatives

  SetCursor(*wxHOURGLASS_CURSOR);
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
    uiState.numComputedNonLinearModes = oldNumNonLinearModes;
    return 3;
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

  // === compute cubic polynomials

  uiState.numComputationThreads = numComputationThreads;

  StVKReducedInternalForces * newCubicPolynomials = NULL;

  SetCursor(*wxHOURGLASS_CURSOR);
  CubicPolynomialsWorker(&code, &newCubicPolynomials);
  SetCursor(*wxSTANDARD_CURSOR);

  if (code != 0)
  {
    this->errMsg( _T("Cubic polynomials computation failed"),
           _T("Cubic polynomials computation failed.") );
    free(newCubicPolynomials);
    return 4;
  }
  else
  {
    delete(precomputationState.cubicPolynomials);
    precomputationState.cubicPolynomials = newCubicPolynomials;
    precomputationState.cubicPolynomialsAvailable = true;

    UpdateMenus();

    Refresh();
  }

  // === save the resulting nonlinear modes and cub files
  if (uiState.enableBatchOutput)
  {
    // save nonlinear modes
    wxString nonLinearModesFilename = projectName + _T(".U");
    const char * nonLinearModesFilenameC = (const char*)nonLinearModesFilename.mb_str();
    int code = WriteMatrixToDisk((char*)nonLinearModesFilenameC,
      3 * precomputationState.nonLinearModalMatrix->Getn(),
      precomputationState.nonLinearModalMatrix->Getr(),
      precomputationState.nonLinearModalMatrix->GetMatrix());
    code = code + 1;

    // save cubic polynomials
    wxString cubicPolynomialsFilename = projectName + _T(".cub");
    const char * cubicPolynomialsFilenameC = (const char*)cubicPolynomialsFilename.mb_str();
    code = precomputationState.cubicPolynomials->Save((char*)cubicPolynomialsFilenameC);
  }

  return 0;
}

void MyFrame::OnBatchPreprocessSingleSimulation(wxCommandEvent & event)
{
  // === create the dialog

  wxDialog * dlg = new wxDialog(this, -1, _T("Batch pre-process single model"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );

  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  // === num linear modes control

  // set suggested number of linear modes
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

  wxStaticText * numLinearModesText = new wxStaticText(dlg, -1,
       _T("Number of linear modes: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numLinearModesControl = new wxTextCtrl(dlg, -1,
      wxString(numLinearModesStringC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numLinearModesSizer = new wxBoxSizer(wxHORIZONTAL);
  numLinearModesSizer->Add(numLinearModesText, 0, wxCENTER);
  numLinearModesSizer->Add(numLinearModesControl, 0, wxCENTER);
  dlgSizer->Add(numLinearModesSizer, 0, wxALL | wxALIGN_LEFT, 4);

  // === num nonlinear modes control

  // set suggested number of nonlinear modes
  char numNonLinearModesStringC[256];
  sprintf(numNonLinearModesStringC, "%d", uiState.numComputedNonLinearModes);

  wxStaticText * numNonLinearModesText = new wxStaticText(dlg, -1,
       _T("Size of simulation basis: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numNonLinearModesControl = new wxTextCtrl(dlg, -1,
      wxString(numNonLinearModesStringC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numNonLinearModesSizer = new wxBoxSizer(wxHORIZONTAL);
  numNonLinearModesSizer->Add(numNonLinearModesText, 0, wxCENTER);
  numNonLinearModesSizer->Add(numNonLinearModesControl, 0, wxCENTER);
  dlgSizer->Add(numNonLinearModesSizer, 0, wxALL | wxALIGN_LEFT, 4);

  // === num computation threads control

  // set the number of computation threads
  char numComputationThreadsC[96];
  sprintf(numComputationThreadsC, "%d", uiState.numComputationThreads);

  wxStaticText * numComputationThreadsText = new wxStaticText(dlg, -1,
       _T("Number of computation threads: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numComputationThreadsControl = new wxTextCtrl(dlg, -1,
      wxString(numComputationThreadsC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numComputationThreadsSizer = new wxBoxSizer(wxHORIZONTAL);
  numComputationThreadsSizer->Add(numComputationThreadsText, 0, wxCENTER);
  numComputationThreadsSizer->Add(numComputationThreadsControl, 0, wxCENTER);
  dlgSizer->Add(numComputationThreadsSizer, 0, wxALL | wxALIGN_LEFT, 4);

  // === output project name control

  wxCheckBox * enableBatchOutput = new wxCheckBox(dlg, -1, _T("Output modal matrix and polynomials to disk files"));
  enableBatchOutput->SetValue(uiState.enableBatchOutput);
  dlgSizer->Add(enableBatchOutput, 0, wxALL | wxALIGN_LEFT, 4);

  wxStaticText * projectNameText = new wxStaticText(dlg, -1,
       _T("Output prefix: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * projectNameControl = new wxTextCtrl(dlg, -1, uiState.currentWorkingDirectory + _T("/") + uiState.projectName, wxDefaultPosition, wxSize(400, -1));

  wxBoxSizer * projectNameSizer = new wxBoxSizer(wxHORIZONTAL);
  projectNameSizer->Add(projectNameText, 0, wxALIGN_CENTER);
  projectNameSizer->Add(projectNameControl, 0, wxEXPAND | wxCENTER);
  dlgSizer->Add(projectNameSizer, 0, wxLEFT | wxRIGHT | wxEXPAND | wxALIGN_LEFT, 4);

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxALL | wxALIGN_CENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  uiState.enableBatchOutput = enableBatchOutput->GetValue();
  wxString numLinearModes_valueString = numLinearModesControl->GetValue();
  wxString numNonLinearModes_valueString = numNonLinearModesControl->GetValue();
  wxString numComputationThreads_valueString = numComputationThreadsControl->GetValue();
  wxString projectName = projectNameControl->GetValue();

  delete(dlg);

  // check if user input is legal

  bool goodInput = true;

  long numLinearModes_value;
  goodInput = goodInput & numLinearModes_valueString.ToLong(&numLinearModes_value);
  long numNonLinearModes_value;
  goodInput = goodInput & numNonLinearModes_valueString.ToLong(&numNonLinearModes_value);
  long numComputationThreads_value;
  goodInput = goodInput & numComputationThreads_valueString.ToLong(&numComputationThreads_value);

  if (goodInput)
  {
    if ((numLinearModes_value <= 0) || (numLinearModes_value > 16384))
      goodInput = false;

    if ((numNonLinearModes_value <= 0) || (numNonLinearModes_value > 16384))
      goodInput = false;

    if ((numComputationThreads_value <= 0) || (numComputationThreads_value > 65536))
      goodInput = false;
  }

  if (goodInput)
  {
    int code = PreprocessSingleSimulation(projectName, numLinearModes_value, numNonLinearModes_value, numComputationThreads_value);
    if (code != 0)
      printf("Error: preprocess failed, exit code %d.\n", code);
  }
  else
  {
    this->errMsg( _T("Invalid input"),  _T("Invalid input") );
  }
}

void MyFrame::PreprocessManySimulationsDeallocateHelper(int numModels, char ** meshFiles, char ** bouFiles, int * numLinearModes, int * numNonLinearModes)
{
  free(numNonLinearModes);
  free(numLinearModes);
  for(int i=0; i<numModels; i++)
  {
    free(meshFiles[i]);
    free(bouFiles[i]);
  }

  free(meshFiles);
  free(bouFiles);
}

class PreprocessManySimulationsDialog : public wxDialog
{
public:
  PreprocessManySimulationsDialog(wxWindow * parent, MyFrame * myFrame);

  wxString GetScriptFilename() { return scriptFilenameTextCtrl->GetValue(); }
  bool GetNumComputationThreads(long * numComputationThreads) { return numComputationThreadsTextCtrl->GetValue().ToLong(numComputationThreads); } 

protected:
  MyFrame * myFrame;

  wxBoxSizer * dlgSizer;

  wxStaticText * scriptFilenameText;
  wxTextCtrl * scriptFilenameTextCtrl;
  wxBitmapButton * scriptFilenameButton;
  wxBoxSizer * scriptFilenameTextSizer;

  wxStaticText * numComputationThreadsText;
  wxTextCtrl * numComputationThreadsTextCtrl;
  wxBoxSizer * numComputationThreadsSizer;

  void OnScriptFilenameButton(wxCommandEvent& event);

  DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(PreprocessManySimulationsDialog, wxDialog)
    EVT_BUTTON(ID_ScriptFilenameButton, PreprocessManySimulationsDialog::OnScriptFilenameButton)
END_EVENT_TABLE()

PreprocessManySimulationsDialog::PreprocessManySimulationsDialog(wxWindow * parent, MyFrame * myFrame_):
  wxDialog(parent, -1, _T(""), wxDefaultPosition, wxDefaultSize), myFrame(myFrame_)
{
  wxString title = _T("Batch pre-process many models");
  SetTitle(title);

  dlgSizer = new wxBoxSizer(wxVERTICAL);

  // === script filename control

  scriptFilenameText = new wxStaticText(this, -1,
       _T("Input script: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxString defaultScriptFilename = myFrame->uiState.manySimulationsscriptFilename;
  if (defaultScriptFilename == _T("none"))
    defaultScriptFilename = myFrame->uiState.currentWorkingDirectory + _T("/");

  scriptFilenameTextCtrl = new wxTextCtrl(this, -1, defaultScriptFilename, wxDefaultPosition, wxSize(400, -1));

  wxImage openFileImage(folder_open_xpm_2);
  scriptFilenameButton = new wxBitmapButton(
    this, ID_ScriptFilenameButton, wxBitmap(openFileImage));

  scriptFilenameTextSizer = new wxBoxSizer(wxHORIZONTAL);
  scriptFilenameTextSizer->Add(scriptFilenameText, 0, wxALIGN_CENTER);
  scriptFilenameTextSizer->Add(scriptFilenameTextCtrl, 0, wxEXPAND | wxCENTER);
  scriptFilenameTextSizer->Add(scriptFilenameButton, 0, wxCENTER);
  dlgSizer->Add(scriptFilenameTextSizer, 0, wxLEFT | wxRIGHT | wxEXPAND | wxALIGN_LEFT, 4);

  // === num computation threads control

  char numComputationThreadsC[96];
  sprintf(numComputationThreadsC, "%d", myFrame->uiState.numComputationThreads);

  numComputationThreadsText = new wxStaticText(this, -1,
       _T("Number of computation threads: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  numComputationThreadsTextCtrl = new wxTextCtrl(this, -1,
      wxString(numComputationThreadsC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  numComputationThreadsSizer = new wxBoxSizer(wxHORIZONTAL);
  numComputationThreadsSizer->Add(numComputationThreadsText, 0, wxCENTER);
  numComputationThreadsSizer->Add(numComputationThreadsTextCtrl, 0, wxCENTER);
  dlgSizer->Add(numComputationThreadsSizer, 0, wxALL | wxALIGN_LEFT, 4);

  wxSizer * buttonSizer = CreateButtonSizer(wxOK | wxCANCEL);
  dlgSizer->Add(buttonSizer, 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(this);
}


void PreprocessManySimulationsDialog::OnScriptFilenameButton(wxCommandEvent& event)
{
  wxString title = _T("Select script filename");
  wxString fileTypes = _T("All files(*)|*");

  wxFileDialog *dlg = new wxFileDialog(this, title, myFrame->uiState.currentWorkingDirectory,
                _T(""), fileTypes,
                wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString filename(dlg->GetPath());
    myFrame->SaveCurrentWorkingDirectory(filename);
    if( !filename.empty() )
      scriptFilenameTextCtrl->SetValue(filename);
  }
}

void MyFrame::OnBatchPreprocessManySimulations(wxCommandEvent & event)
{
/*
  // === create the dialog

  wxDialog * dlg = new wxDialog(this, -1, _T("Batch pre-process many models"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );

  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  // === script filename control

  wxStaticText * scriptFilenameText = new wxStaticText(dlg, -1,
       _T("Input script: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * scriptFilenameControl = new wxTextCtrl(dlg, -1, uiState.currentWorkingDirectory + _T("/"), wxDefaultPosition, wxSize(400, -1));

  wxImage openFileImage(folder_open_xpm_2);
  wxBitmapButton * scriptFilenameButton = new wxBitmapButton(
    this, ID_ScriptFilenameButton, wxBitmap(openFileImage));

  wxBoxSizer * projectNameSizer = new wxBoxSizer(wxHORIZONTAL);
  projectNameSizer->Add(scriptFilenameText, 0, wxALIGN_CENTER);
  projectNameSizer->Add(scriptFilenameControl, 0, wxEXPAND | wxCENTER);
  projectNameSizer->Add(scriptFilenameButton, 0, wxCENTER);
  dlgSizer->Add(projectNameSizer, 0, wxLEFT | wxRIGHT | wxEXPAND | wxALIGN_LEFT, 4);

  // === num computation threads control

  char numComputationThreadsC[96];
  sprintf(numComputationThreadsC, "%d", uiState.numComputationThreads);

  wxStaticText * numComputationThreadsText = new wxStaticText(dlg, -1,
       _T("Number of computation threads: "), wxDefaultPosition, wxDefaultSize,
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * numComputationThreadsControl = new wxTextCtrl(dlg, -1,
      wxString(numComputationThreadsC, wxConvUTF8), wxDefaultPosition, wxSize(100,-1));

  wxBoxSizer * numComputationThreadsSizer = new wxBoxSizer(wxHORIZONTAL);
  numComputationThreadsSizer->Add(numComputationThreadsText, 0, wxCENTER);
  numComputationThreadsSizer->Add(numComputationThreadsControl, 0, wxCENTER);
  dlgSizer->Add(numComputationThreadsSizer, 0, wxALL | wxALIGN_LEFT, 4);

  // === OK, CANCEL buttons

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxALL | wxALIGN_CENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString numComputationThreads_valueString = numComputationThreadsControl->GetValue();
  wxString scriptFilename = scriptFilenameControl->GetValue();

  delete(dlg);

  long numComputationThreads_value;
  bool goodInput = numComputationThreads_valueString.ToLong(&numComputationThreads_value);
*/

  PreprocessManySimulationsDialog * dlg = new PreprocessManySimulationsDialog(this, this);
  
  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  printf("Here...\n"); fflush(NULL);

  long numComputationThreads_value;
  bool goodInput = dlg->GetNumComputationThreads(&numComputationThreads_value);
  wxString scriptFilename = dlg->GetScriptFilename();
  delete(dlg);

  // check if user input is legal
  if (!goodInput)
  {
    this->errMsg( _T("Incorrect number of threads."), _T("Incorrect number of threads."));
    return;
  }

  uiState.numComputationThreads = numComputationThreads_value;
  uiState.manySimulationsscriptFilename = scriptFilename;

  printf("Num computation threads: %ld\n", numComputationThreads_value); fflush(NULL);
  printf("Opening script file %s...\n", (const char*) scriptFilename.mb_str()); fflush(NULL);

  FILE * fin = fopen((const char*)scriptFilename.mb_str(), "r");
  if (!fin)
  {
    this->errMsg( _T("Could not open the script file."), _T("Could not open the script file."));
    return;
  }

  // read the number of models
  int numModels;
  if ((fscanf(fin, "%d\n", &numModels) < 0) || (numModels <= 0))
  {
    this->errMsg( _T("Incorrect number of models in the script file."), _T("Incorrect number of models in the script file."));
    fclose(fin);
    return;
  }

  printf("Detected %d models.\n", numModels);

  // read directory name
  char directoryName[4096];
  if (fgets(directoryName, 4096, fin) == NULL)
  {
    this->errMsg( _T("Incorrect directory name in the script file."), _T("Incorrect directory name in the script file."));
    fclose(fin);
    return;
  }

  if (directoryName[strlen(directoryName)-1] == '\n')
    directoryName[strlen(directoryName)-1] = 0;

  printf("Directory name is %s .\n", directoryName);

  char ** meshFiles = (char**) malloc (sizeof(char*) * numModels);
  char ** bouFiles = (char**) malloc (sizeof(char*) * numModels);
  for(int i=0; i<numModels; i++)
  {
    meshFiles[i] = (char*) malloc (sizeof(char) * 256);
    bouFiles[i] = (char*) malloc (sizeof(char) * 256);
  }
  int * numLinearModes = (int*) malloc (sizeof(int) * numModels);
  int * numNonLinearModes = (int*) malloc (sizeof(int) * numModels);
  
  for(int i=0; i<numModels; i++)
  {
    if (fscanf(fin, "%s %s %d %d\n", meshFiles[i], bouFiles[i], &numLinearModes[i], &numNonLinearModes[i]) != 4)
    {
      wxString errM = _T("Incorrect line in the script file, model #");
      char numStr[128] = "";
      sprintf(numStr,"%i .",i);
      errM += wxString(numStr, wxConvUTF8);
      this->errMsg(errM, errM);
      fclose(fin);
      PreprocessManySimulationsDeallocateHelper(numModels, meshFiles, bouFiles, numLinearModes, numNonLinearModes);
      return;
    }
    else
    {
      printf("Model #%d: %s %s %d %d\n", i, meshFiles[i], bouFiles[i], numLinearModes[i], numNonLinearModes[i]);
    }
  }
  fclose(fin);

  // check that the number of linear and nonlinear modes are meaningful
  for(int i=0; i<numModels; i++)
  {
    if ((numLinearModes[i] <= 0) || (numLinearModes[i] > 16384))
    {
      wxString errM = _T("Incorrect number of linear modes for model #");
      char numStr[128] = "";
      sprintf(numStr,"%i .",i);
      errM += wxString(numStr, wxConvUTF8);
      this->errMsg(errM, errM);
      PreprocessManySimulationsDeallocateHelper(numModels, meshFiles, bouFiles, numLinearModes, numNonLinearModes);
      return;
    }

    if ((numNonLinearModes[i] <= 0) || (numNonLinearModes[i] > 16384))
    {
      wxString errM = _T("Incorrect number of nonlinear modes for model #");
      char numStr[128] = "";
      sprintf(numStr,"%i .",i);
      errM += wxString(numStr, wxConvUTF8);
      this->errMsg(errM, errM);
      PreprocessManySimulationsDeallocateHelper(numModels, meshFiles, bouFiles, numLinearModes, numNonLinearModes);
      return;
    }
  }

  // do the pre-process
  for(int i=0; i<numModels; i++)
  {
    printf("\n\n");
    printf("*****************************************\n");
    printf("**** Pre-processing model #%d/%d: %s %s numLin=%d numNonLin=%d...\n", i, numModels, meshFiles[i], bouFiles[i], numLinearModes[i], numNonLinearModes[i]);
    printf("*****************************************\n\n");

    wxString meshFilewx = wxString(meshFiles[i], wxConvUTF8);
    wxString localProjectName = meshFilewx.BeforeLast('.');
    wxString projectName = wxString(directoryName,  wxConvUTF8) + _T("/") + localProjectName;

    // load the simulation mesh
    wxString meshFilenameComplete = wxString(directoryName,  wxConvUTF8) + _T("/") + meshFilewx;
    int meshCode = LoadSimulationMesh(meshFilenameComplete);
    if (meshCode != 0)
    {
      char s[96];
      sprintf(s, "Error: preprocess %d failed, could not load mesh %s, exit code %d.\n", i, (const char*)meshFilenameComplete.mb_str(), meshCode);
      printf("%s\n", s);
      wxString errM(s, wxConvUTF8);
      this->errMsg(errM, errM);
      break;
    }

    // load the bou file
    wxString bouFilewx = wxString(bouFiles[i], wxConvUTF8);
    wxString bouFilenameComplete = wxString(directoryName,  wxConvUTF8) + _T("/") + bouFilewx;
    int bouCode = LoadFixedVertices(bouFilenameComplete);
    if (bouCode != 0)
    {
      char s[96];
      sprintf(s, "Error: preprocess %d failed, could not load fixed vertices %s, exit code %d.\n", i, (const char*)bouFilenameComplete.mb_str(), bouCode);
      printf("%s\n", s);
      wxString errM(s, wxConvUTF8);
      this->errMsg(errM, errM);
      break;
    }

    int code = PreprocessSingleSimulation(projectName, numLinearModes[i], numNonLinearModes[i], numComputationThreads_value);
    if (code != 0)
    {
      char s[96];
      sprintf(s, "Error: preprocess %d failed, exit code %d.\n", i, code);
      printf("%s\n", s);
      wxString errM(s, wxConvUTF8);
      this->errMsg(errM, errM);
      break;
    }
  }

  PreprocessManySimulationsDeallocateHelper(numModels, meshFiles, bouFiles, numLinearModes, numNonLinearModes);

  printf("\nPreprocessing completed (%d models).\n", numModels);
}

void MyFrame::OnPrepareRuntimeSimulation(wxCommandEvent& event)
{
  int createdRenderingMeshFromSimulationMesh = 0;
  if (!precomputationState.renderingMeshAvailable)
  {
    wxMessageDialog * dlg = new wxMessageDialog(this, 
      _T(
         "Click OK to use the surface mesh of the simulation mesh as the rendering mesh.\n\n"
         L"Click CANCEL to load your own rendering mesh (via the \"Mesh.Load Triangle Mesh\" menu option)."
        ),
      _T("Rendering mesh has not been specified"),
      wxOK | wxCANCEL);

    if (dlg->ShowModal() != wxID_OK)
    {
      delete(dlg);
      return;
    }

    delete(dlg);

    // create rendering mesh from the simulation mesh
    CreateRenderingMeshFromSimulationMesh();
    createdRenderingMeshFromSimulationMesh = 1;
  }

  // get the project name, and the scale to 1 Hz flag
  wxDialog * dlg = new wxDialog(this, -1, _T("Prepare real-time simulation"),
    wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("dialogBox") );
  
  wxBoxSizer * dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxStaticText * projectNameText = new wxStaticText(dlg, -1, 
       _T("Project name: "), wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER, _T( "staticText"));

  wxTextCtrl * projectNameControl = new wxTextCtrl(dlg, -1, uiState.projectName, wxDefaultPosition, wxSize(200, -1));

  wxBoxSizer * projectNameSizer = new wxBoxSizer(wxHORIZONTAL);
  projectNameSizer->Add(projectNameText, 0, wxCENTER);
  projectNameSizer->Add(projectNameControl, 0, wxCENTER);
  dlgSizer->Add(projectNameSizer, 0, wxALL | wxCENTER, 8);

  wxCheckBox * scaleTo1HzCheckBox = NULL;
  if (precomputationState.linearModesAvailable)
  {
    scaleTo1HzCheckBox = new wxCheckBox(dlg, 
      -1, _T("Scale frequency spectrum to 1 Hz"));
    scaleTo1HzCheckBox->SetValue(uiState.scaleRealTimeTo1HzCheckBox);
    dlgSizer->Add(scaleTo1HzCheckBox, 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);
  }

  dlgSizer->Add(dlg->CreateButtonSizer(wxOK | wxCANCEL), 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  dlg->SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(dlg);

  if (dlg->ShowModal() != wxID_OK)
  {
    if (createdRenderingMeshFromSimulationMesh)
      DeleteRenderingMesh();
    delete(dlg);
    return;
  }

  uiState.projectName = projectNameControl->GetValue();
  if ((precomputationState.linearModesAvailable) && (precomputationState.frequenciesAvailable))
    uiState.scaleRealTimeTo1HzCheckBox = scaleTo1HzCheckBox->GetValue();
  else
    uiState.scaleRealTimeTo1HzCheckBox = false; 
  delete(dlg);

  // create subdirectory
  const char * projectNameC = (const char*)uiState.projectName.mb_str();
  char s[4096];

  printf("\nPreparing a real-time simulation:\n");
  printf("Creating subdirectory: realtime/%s\n", projectNameC);

  // first, create the top level "realtime" directory
  sprintf(s,"realtime");
  if (!wxDirExists(s))
    wxMkdir(_T("realtime"));

  // create the project subdirectory
  sprintf(s,"realtime/%s",projectNameC);
  if (!wxDirExists(s))
    wxMkdir(_T("realtime/") + uiState.projectName);

  // copy rendering mesh file, with mtl
  printf("Writing rendering mesh and rendering materials file.\n");
  sprintf(s,"realtime/%s/%s.obj", projectNameC, projectNameC);
  precomputationState.renderingMesh->save(string(s), 1);
  
  // generate .URendering.float
  if (!precomputationState.interpolationDataAvailable)
    BuildInterpolant();

  printf("Interpolating and saving URendering.float matrix.\n");
  double * inputMatrix = precomputationState.nonLinearModalMatrix->GetMatrix();
  int nTarget = (int)(precomputationState.renderingMesh->getNumVertices());
  int dataSize = 3 * nTarget * precomputationState.nonLinearModalMatrix->Getr();
  double * outputMatrix = (double*) malloc (sizeof(double) * dataSize);
  for(int i=0; i<precomputationState.nonLinearModalMatrix->Getr(); i++)
  {
    precomputationState.simulationMesh->interpolate(
         &inputMatrix[ELT(3*precomputationState.simulationMesh->getNumVertices(),0,i)],
         &outputMatrix[ELT(3*nTarget,0,i)],
         nTarget, precomputationState.simulationMesh->getNumElementVertices(),
         precomputationState.interpolationData_vertices,
         precomputationState.interpolationData_weights);
  }

  float * outputMatrixFloat = (float*) outputMatrix;

  for(int i=0; i< dataSize; i++)
    outputMatrixFloat[i] = (float)(outputMatrix[i]);
  
  // save file to disk
  sprintf(s,"realtime/%s/%s.URendering.float", projectNameC, projectNameC);
  int code = WriteMatrixToDisk(s,
       3 * nTarget, 
       precomputationState.nonLinearModalMatrix->Getr(), 
       outputMatrixFloat);

  if (code != 0)
  {
    printf("Error: failed to write to file: realtime/%s/%s.URendering.float .\n", projectNameC, projectNameC);
  }

  free(outputMatrix); 

  // copy .cub
  printf("Saving cub file.\n");
  sprintf(s,"realtime/%s/%s.cub", projectNameC, projectNameC);
  precomputationState.cubicPolynomials->Save(s);

  // create .lighting
  printf("Creating lighting file.\n");
  sprintf(s,"realtime/%s/%s.lighting", projectNameC, projectNameC);
  FILE * fout;
  OpenFile_(s, &fout, "w");
  fprintf(fout, "%s", lightingString);
  fclose(fout);

  // create .config
  printf("Creating config file.\n");
  sprintf(s,"realtime/%s/%s.config", projectNameC, projectNameC);
  OpenFile_(s, &fout, "w");

  fprintf(fout, "*deformableObjectFilename\n");
  fprintf(fout, "realtime/%s/%s.obj\n\n", projectNameC, projectNameC);

  fprintf(fout, "*modesFilename\n");
  fprintf(fout, "realtime/%s/%s.URendering.float\n\n", projectNameC, projectNameC);

  fprintf(fout, "*cubicPolynomialFilename\n");
  fprintf(fout, "realtime/%s/%s.cub\n\n", projectNameC, projectNameC);

  // set compliance automatically
  SparseMatrix * massMatrix;
  GenerateMassMatrix::computeMassMatrix(precomputationState.simulationMesh, &massMatrix, true);
  double totalMass = massMatrix->SumEntries();
  // compute stiffness matrix
  int r = (precomputationState.nonLinearModalMatrix)->Getr();
  int n3 = 3 * (precomputationState.nonLinearModalMatrix)->Getn();
  double * K = (double*) malloc (sizeof(double) * r * r);
  double * zero = (double*) calloc (r, sizeof(double));

  StVKReducedStiffnessMatrix stVKReducedStiffnessMatrix(precomputationState.cubicPolynomials);
  stVKReducedStiffnessMatrix.Evaluate(zero, K);
  if (uiState.scaleRealTimeTo1HzCheckBox)
  {
    double factor = (precomputationState.frequencies)[precomputationState.numRigidModes];
    factor = 1.0 / (factor * factor);
    for(int i=0; i<r*r; i++)
      K[i] *= factor;
  }
  Matrix<double> KM(r, r, K, false, false);
  Matrix<double> KInverseM = Inverse(KM);
  double * KInvUT = (double*) malloc (sizeof(double) * n3 * r);
  (precomputationState.nonLinearModalMatrix)->AssembleMatrix(r, KInverseM.GetData(), KInvUT);
  InPlaceTransposeMatrix(n3, r, KInvUT);
  double gamma2 = 0;
  // find column of KInvUT with maximum 2-norm
  for(int i=0; i<n3; i++) // over all columns
  {
    double norm2 = 0.0;
    for(int j=0; j<r; j++)
      norm2 += KInvUT[ELT(r,j,i)] * KInvUT[ELT(r,j,i)];
    if (norm2 > gamma2)
      gamma2 = norm2;
  }

  double gamma = sqrt(gamma2);
  double radius = precomputationState.simulationMeshGeometricParameters.radius;
  double beta = 0.25;
  double numPixels = 100.0;
  double compliance = sqrt(totalMass) * beta * radius / (gamma * numPixels);

  fprintf(fout, "*deformableObjectCompliance\n");
  fprintf(fout, "%G\n\n", compliance);

  free(KInvUT);
  free(zero);
  free(K);
  delete(massMatrix);

  fprintf(fout, "*frequencyScaling\n");
  if (uiState.scaleRealTimeTo1HzCheckBox)
    fprintf(fout, "%G\n\n", 1.0 / (precomputationState.frequencies)[precomputationState.numRigidModes]);
  else
    fprintf(fout, "1.0\n\n");

  fprintf(fout, "*dampingMassCoef\n");
  fprintf(fout, "0.0\n\n");

  fprintf(fout, "*dampingStiffnessCoef\n");
  fprintf(fout, "0.003\n\n");

  fprintf(fout, "*substepsPerTimeStep\n");
  fprintf(fout, "5\n\n");

  fprintf(fout, "*renderWireframe\n");
  fprintf(fout, "1\n\n");

  fprintf(fout, "*cameraRadius\n");
  fprintf(fout, "%f\n\n", myGLCanvas->GetCamera()->GetRadius());

  fprintf(fout, "*cameraLattitude\n");
  fprintf(fout, "%f\n\n", myGLCanvas->GetCamera()->GetTheta() * 180 / M_PI);

  fprintf(fout, "*cameraLongitude\n");
  fprintf(fout, "%f\n\n", myGLCanvas->GetCamera()->GetPhi() * 180 / M_PI);

  double focusPosition[3];
  myGLCanvas->GetCamera()->GetFocusPosition(focusPosition);

  fprintf(fout, "*focusPositionX\n");
  fprintf(fout, "%f\n\n", focusPosition[0]);

  fprintf(fout, "*focusPositionY\n");
  fprintf(fout, "%f\n\n", focusPosition[1]);

  fprintf(fout, "*focusPositionZ\n");
  fprintf(fout, "%f\n\n", focusPosition[2]);

  fprintf(fout, "*backgroundColor\n");
  fprintf(fout, "220 220 220\n\n");

  fprintf(fout, "*lightingConfigFilename\n");
  fprintf(fout, "realtime/%s/%s.lighting\n\n", projectNameC, projectNameC);

  fprintf(fout, "*renderOnGPU\n");
  fprintf(fout, "1\n\n");

  fprintf(fout, "*windowWidth\n");
  fprintf(fout, "800\n\n");

  fprintf(fout, "*windowHeight\n");
  fprintf(fout, "800\n\n");

  fclose(fout);

  // create .bat file
  sprintf(s,"%s.bat", projectNameC);
  printf("Creating batch file %s.\n", s);
  OpenFile_(s, &fout, "w");

  fprintf(fout, "../reducedDynamicSolver-rt/reducedDynamicSolver-rt realtime/%s/%s.config", projectNameC, projectNameC);
  fclose(fout);

  #ifdef __WXMSW__
  #else
    sprintf(s,"chmod u+x %s.bat", projectNameC);
    int systemCode = system(s);
    if (systemCode != 0)
      printf("Error: exit code from system call: %d\n", systemCode);
  #endif

  precomputationState.runtimeSimReady = true;

  if (createdRenderingMeshFromSimulationMesh)
    DeleteRenderingMesh();

  UpdateMenus();

  Refresh();
}

void MyFrame::OnLaunchRuntimeSimulation(wxCommandEvent& event)
{
  // launch external runtime simulation
  const char * projectNameC = (const char*)uiState.projectName.mb_str();
  char s[4096];
  #ifdef __WXMSW__
    sprintf(s,"%s.bat", projectNameC);
  #else
    sprintf(s,"./%s.bat", projectNameC);
  #endif
  int systemCode = system(s);
  if (systemCode != 0)
    printf("Error: exit code from system call: %d\n", systemCode);
}

void MyFrame::OnTuneRuntimeSimulation(wxCommandEvent& event)
{
  wxString message = 
	  _T( 
  "You can tune the real-time simulation from within the real-time GUI (activated by choosing \"Launch real-time simulation\").\n\n"
  L"There are essentially two parameters to tune:\n"
  L"  1. You can set how stiff the object feels by adjusting the \"strength\" of your mouse forces via the \"Deformable object compliance\" parameter.\n"
  L"  2. You can make the object oscillate faster/slower by altering the \"Base frequency\" parameter.\n"
   );
        
  wxMessageBox(message,
    _T("Tuning the real-time simulation"), wxOK | wxICON_INFORMATION, this);
}

