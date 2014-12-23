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

// matrix conversion routines

#include "largeModalDeformationFactory.h"
#include "matrixIO.h"

/* XPM */
static char *folder_open_xpm[] = {
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
"                ",
".....           ",
".XXXo.          ",
".XXXXO........  ",
".XXXXXXXXXXXX.  ",
".XXXXXXXXXXXX.  ",
".X++++++++++++++",
".X+@#$%&*=-;:>,+",
".<.1@#$%2*=-;:23",
"..X41@#$%2*=-;3 ",
"..X561@#$%2*=-3 ",
".78X561@#$%2*%3 ",
"90<8X561@#$%23  ",
"q++++++++++++w  ",
"                "
};

class MatrixConversionDialog : public wxDialog
{
public:
  MatrixConversionDialog(wxWindow * parent, int dialogType, wxString & defaultWorkingDirectory);

  wxString GetInputMatrix() { return inputMatrixTextCtrl->GetValue(); }
  wxString GetOutputMatrix() { return outputMatrixTextCtrl->GetValue(); }
  wxString GetSelectedFilename() { return selectedFilename; }

protected:
  int dialogType; // 0 = text to binary, 1 = binary to text

  wxBoxSizer * dlgSizer;

  wxStaticText * inputMatrixText;
  wxTextCtrl * inputMatrixTextCtrl;
  wxBitmapButton * inputMatrixFilenameButton;
  wxBoxSizer * inputMatrixTextSizer;

  wxStaticText * outputMatrixText;
  wxTextCtrl * outputMatrixTextCtrl;
  wxBitmapButton * outputMatrixFilenameButton;
  wxBoxSizer * outputMatrixTextSizer;

  wxString defaultWorkingDirectory;
  wxString selectedFilename;

  void OnInputMatrixFilenameButton(wxCommandEvent& event);
  void OnOutputMatrixFilenameButton(wxCommandEvent& event);

  void SelectMatrixFilename(int dialogType, int matrixType, wxTextCtrl * targetCtrl);
   
  DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(MatrixConversionDialog, wxDialog)
    EVT_BUTTON(ID_inputMatrixFilenameButton, MatrixConversionDialog::OnInputMatrixFilenameButton)
    EVT_BUTTON(ID_outputMatrixFilenameButton, MatrixConversionDialog::OnOutputMatrixFilenameButton)
END_EVENT_TABLE()

MatrixConversionDialog::MatrixConversionDialog(wxWindow * parent, int dialogType, wxString & defaultWorkingDirectory_):
  wxDialog(parent, -1, _T(""), wxDefaultPosition, wxDefaultSize), defaultWorkingDirectory(defaultWorkingDirectory_), selectedFilename(_T(""))
{
  this->dialogType = dialogType;
  wxString title;
  if (dialogType == 0)
    title = _T("Import matrix from text file");
  else
    title = _T("Export matrix to text file");
  SetTitle(title);

  dlgSizer = new wxBoxSizer(wxVERTICAL);

  wxString inputMatrixTextString = _T("Input ");
  if (dialogType == 0) 
    inputMatrixTextString.append(_T("text"));
  else
    inputMatrixTextString.append(_T("binary"));
  inputMatrixTextString.append(_T(" matrix file:"));

  inputMatrixText = new wxStaticText(this, -1, 
       inputMatrixTextString, wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER);

  inputMatrixTextCtrl = new wxTextCtrl(this, -1, 
      wxEmptyString, wxDefaultPosition, wxSize(200,-1), wxTE_LEFT);

  wxImage openFileImage(folder_open_xpm);

  inputMatrixFilenameButton = new wxBitmapButton(
    this, ID_inputMatrixFilenameButton, wxBitmap(openFileImage));

  inputMatrixTextSizer = new wxBoxSizer(wxHORIZONTAL);
  inputMatrixTextSizer->Add(inputMatrixText, 0, wxCENTER);
  inputMatrixTextSizer->Add(inputMatrixTextCtrl, 0, wxCENTER);
  inputMatrixTextSizer->Add(inputMatrixFilenameButton, 0, wxCENTER);
  dlgSizer->Add(inputMatrixTextSizer, 0, wxALL | wxCENTER, 8);

  wxString outputMatrixTextString = _T("Output ");
  if (dialogType == 0) 
    outputMatrixTextString.append(_T("binary"));
  else
    outputMatrixTextString.append(_T("text"));
  outputMatrixTextString.append(_T(" matrix file:"));

  outputMatrixText = new wxStaticText(this, -1, 
       outputMatrixTextString, wxDefaultPosition, wxDefaultSize, 
       wxALIGN_CENTER);

  outputMatrixTextCtrl = new wxTextCtrl(this, -1, 
      wxEmptyString, wxDefaultPosition, wxSize(200,-1), wxTE_LEFT);

  outputMatrixFilenameButton = new wxBitmapButton(
    this, ID_outputMatrixFilenameButton, wxBitmap(openFileImage));

  outputMatrixTextSizer = new wxBoxSizer(wxHORIZONTAL);
  outputMatrixTextSizer->Add(outputMatrixText, 0, wxCENTER);
  outputMatrixTextSizer->Add(outputMatrixTextCtrl, 0, wxCENTER);
  outputMatrixTextSizer->Add(outputMatrixFilenameButton, 0, wxCENTER);
  dlgSizer->Add(outputMatrixTextSizer, 0, wxALL | wxCENTER, 8);

  wxSizer * buttonSizer = CreateButtonSizer(wxOK | wxCANCEL);
  dlgSizer->Add(buttonSizer, 0, wxLEFT | wxRIGHT | wxBOTTOM | wxCENTER, 8);

  SetSizer(dlgSizer);
  dlgSizer->SetSizeHints(this);
}

// matrixType:
//   0 = text
//   1 = binary
void MatrixConversionDialog::SelectMatrixFilename(int dialogType, int matrixType, wxTextCtrl * targetCtrl)
{
  wxString title;
  wxString fileTypes;
  if (matrixType == 0)
  {
    if (dialogType == 0)
      title = _T("Select input text matrix");
    else
      title = _T("Select output text matrix");
    fileTypes = _T("Text files(*.txt)|*.txt|All files(*.*)|*.*");
  }
  else
  {
    if (dialogType == 0)
      title = _T("Select output binary matrix");
    else
      title = _T("Select input binary matrix");
    fileTypes = _T("Binary files(*.U)|*.U|All files(*.*)|*.*");
  }

  wxFileDialog *dlg = new wxFileDialog(this, title, defaultWorkingDirectory,
		_T(""), fileTypes,
		wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString filename(dlg->GetPath());
    selectedFilename = filename;
    //SaveCurrentWorkingDirectory(filename);
    if( !filename.empty() )
      targetCtrl->SetValue(filename);
  }
}

void MatrixConversionDialog::OnInputMatrixFilenameButton(wxCommandEvent& event)
{
  if (dialogType == 0)
    SelectMatrixFilename(dialogType, 0, inputMatrixTextCtrl);
  else
    SelectMatrixFilename(dialogType, 1, inputMatrixTextCtrl);
}

void MatrixConversionDialog::OnOutputMatrixFilenameButton(wxCommandEvent& event)
{
  if (dialogType == 0)
    SelectMatrixFilename(dialogType, 1, outputMatrixTextCtrl);
  else
    SelectMatrixFilename(dialogType, 0, outputMatrixTextCtrl);
}

void MyFrame::OnConvertTextMatrixToBinaryMatrix(wxCommandEvent& event)
{
  MatrixConversionDialog * dlg = new MatrixConversionDialog(this, 0, uiState.currentWorkingDirectory);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString selectedFilename = dlg->GetSelectedFilename();
  SaveCurrentWorkingDirectory(selectedFilename);

  wxString inputMatrixFilename =  dlg->GetInputMatrix();
  wxString outputMatrixFilename =  dlg->GetOutputMatrix();

  delete(dlg);

  bool goodInput = (inputMatrixFilename.Length() > 0) && (outputMatrixFilename.Length() > 0);
  if (goodInput)
  {
    int m,n;
    double * U;

    const char * inputFilename = inputMatrixFilename.mb_str();
    int code = ReadMatrixFromDiskTextFile((char*)inputFilename, &m, &n, &U);
    if (code != 0)
    {
      this->errMsg( _T("Loading error"),  
        _T("Unable to load matrix from text file ") + inputMatrixFilename );
      free(U);
      return;
    }

    const char * outputFilename = outputMatrixFilename.mb_str();
    code = WriteMatrixToDisk((char*)outputFilename, m, n, U);
    if (code != 0)
    {
      this->errMsg( _T("Writing error"),  
        _T("Unable to write binary matrix to ") + outputMatrixFilename );
      free(U);
      return;
    }

    free(U);
  }
  else
  {
    printf("Input matrix filename or output matrix filename length was zero.\n");
  }

}

void MyFrame::OnConvertBinaryMatrixToTextMatrix(wxCommandEvent& event)
{
  MatrixConversionDialog * dlg = new MatrixConversionDialog(this, 1, uiState.currentWorkingDirectory);

  if (dlg->ShowModal() != wxID_OK)
  {
    delete(dlg);
    return;
  }

  wxString selectedFilename = dlg->GetSelectedFilename();
  SaveCurrentWorkingDirectory(selectedFilename);

  wxString inputMatrixFilename =  dlg->GetInputMatrix();
  wxString outputMatrixFilename =  dlg->GetOutputMatrix();

  delete(dlg);

  bool goodInput = (inputMatrixFilename.Length() > 0) && (outputMatrixFilename.Length() > 0);
  if (goodInput)
  {
    int m,n;
    double * U;

    const char * inputFilename = inputMatrixFilename.mb_str();
    int code = ReadMatrixFromDisk((char*)inputFilename, &m, &n, &U);
    if (code != 0)
    {
      this->errMsg( _T("Loading error"),  
        _T("Unable to load matrix from binary file ") + inputMatrixFilename );
      free(U);
      return;
    }

    const char * outputFilename = outputMatrixFilename.mb_str();
    code = WriteMatrixToDiskTextFile((char*)outputFilename, m, n, U);
    if (code != 0)
    {
      this->errMsg( _T("Writing error"),  
        _T("Unable to write text matrix to ") + outputMatrixFilename );
      free(U);
      return;
    }

    free(U);
  }
  else
  {
    printf("Input matrix filename or output matrix filename length was zero.\n");
  }
}

