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

// selection of fixed vertices

#include "largeModalDeformationFactory.h"

int MyFrame::LoadFixedVertices(wxString & fixedVerticesFilename)
{
  int * newFixedVertices = NULL;
  int newNumFixedVertices;

  SetCursor(*wxHOURGLASS_CURSOR);
  const char * filename = fixedVerticesFilename.mb_str();
  int code = LoadList::load((char*) filename, &newNumFixedVertices,&newFixedVertices);
  SetCursor(*wxSTANDARD_CURSOR);

  if (code != 0)
  {
    this->errMsg( _T("Loading error"),  
      _T("Unable to load fixed vertices from ") + fixedVerticesFilename );
    free(newFixedVertices);
    return code;
  }

  int fixedVerticedOK = -1;
  int numSimulationVertices = precomputationState.simulationMesh->getNumVertices();
  for(int i=0; i<newNumFixedVertices; i++)
  {
    if ((newFixedVertices[i] < 1) || (newFixedVertices[i] > numSimulationVertices))
    {
      fixedVerticedOK = i;
      break;
    }
  }

  if (fixedVerticedOK != -1)
  {
    char s[4096];
    sprintf(s, "Encountered a fixed vertex that is not a 1-indexed simulation mesh vertex: %d", newFixedVertices[fixedVerticedOK]);
    this->errMsg( _T("Loading error"), wxString(s, wxConvUTF8)); 
    free(newFixedVertices);
    return -1;
  }

  // success
  // convert to 0-indexed
  precomputationState.fixedVertices.clear();
  for(int i=0; i<newNumFixedVertices; i++)
    precomputationState.fixedVertices.insert(newFixedVertices[i] - 1);
  free(newFixedVertices);

  precomputationState.fixedVerticesAvailable = true;
  UpdateMenus();

  Refresh();

  return 0;
}

void MyFrame::OnLoadFixedVertices(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load fixed vertices"),
    uiState.currentWorkingDirectory, _T(""), _T("Comma-separated 1-indexed file(*.bou)|*.bou|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString fixedVerticesFilename( dlg->GetPath() );
    delete(dlg);
    SaveCurrentWorkingDirectory(fixedVerticesFilename);
    if( !fixedVerticesFilename.empty() )
    {
      int code = LoadFixedVertices(fixedVerticesFilename);
      if (code != 0)
      {
        printf("Error: failed to load fixed vertices from file %s , exit code %d.\n", (const char*)fixedVerticesFilename.mb_str(), code);
        return;
      }
    }
  }
  else
  {
    dlg->Destroy();
  }
}

void MyFrame::OnSaveFixedVertices(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Save fixed vertices"),
    uiState.currentWorkingDirectory, _T(""), _T("Fixed Vertices Files(*.bou)|*.bou|All files(*.*)|*.*"),
    wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString vertexFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(vertexFilename);
    if( !vertexFilename.empty() )
    {
      int numFixedVertices = (int)(precomputationState.fixedVertices.size());
      int * fixedVerticesC = (int*) malloc (sizeof(int) * numFixedVertices);
      set<int> :: iterator iter;
      int i=0;
      for(iter = precomputationState.fixedVertices.begin(); iter != precomputationState.fixedVertices.end(); iter++)
      {
        fixedVerticesC[i] = *iter;
        i++;
      }

      //loadList.sort(numFixedVertices, fixedVertices);
      int offset = 1;
      const char * filename = vertexFilename.mb_str();
      int code = LoadList::save((char*)filename, numFixedVertices, fixedVerticesC, offset);
      free(fixedVerticesC);
      if (code != 0)
        this->errMsg( _T("Saving error"), _T("Unable to save fixed vertices to ") +  vertexFilename );
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSelectFixedVertices(wxCommandEvent& event)
{
  ActivateVertexSelection(event.IsChecked());
}

void MyFrame::OnClearFixedVertices(wxCommandEvent& event)
{
  precomputationState.fixedVertices.clear();
  precomputationState.fixedVerticesAvailable = false;
  UpdateMenus();

  Refresh();
}

void MyFrame::OnFixedVerticesInformation(wxCommandEvent& event)
{
  char infoString[4096];

  sprintf(infoString, "  Num fixed vertices: %d\n", (int)precomputationState.fixedVertices.size());  

  wxWindow * parent = NULL;
  WxDialogWrapper<wxMessageDialog> dlg( new wxMessageDialog(parent,
    wxString(infoString, wxConvUTF8), 
    wxString(_T("Fixed Vertices Information"), wxConvUTF8),
      wxOK | wxICON_INFORMATION) );
  dlg.get()->ShowModal();
}

