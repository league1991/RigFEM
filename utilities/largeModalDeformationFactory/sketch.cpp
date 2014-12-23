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

// basis from external deformation data 

#include "matrixIO.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnLoadSketchData(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load external simulation data"), uiState.currentWorkingDirectory, _T(""), _T("External Data Files(*.UData)|*.UData|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString sketchDataFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(sketchDataFilename);
    if( !sketchDataFilename.empty() )
    {
      int newSketchDatasamples;
      double * newSketchData = NULL;

      int n1;
      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = sketchDataFilename.mb_str();
      printf("Loading external simulation data from %s.\n", filename);
      int code = ReadMatrixFromDisk((char*)filename, 
        &n1, &newSketchDatasamples, &newSketchData);
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  
          _T("Unable to load external simulation data from ") + sketchDataFilename );
        dlg->Destroy();
        return;
      }

      if (n1 != 3 * precomputationState.simulationMesh->getNumVertices())
      {
        this->errMsg( _T("Loading error"),  
          _T("The number of vertices in ") + sketchDataFilename + _T(" does not match the simulation mesh."));
        free(newSketchData);
        dlg->Destroy();
        return;
      }

      // success
      delete(precomputationState.sketchDataMatrix);
      precomputationState.sketchDataMatrix = new 
        ModalMatrix(precomputationState.simulationMesh->getNumVertices(), 
          newSketchDatasamples, newSketchData);
      free(newSketchData);

      precomputationState.sketchDataAvailable = true;

      UpdateMenus();

      Refresh();
    }
  }

  dlg->Destroy();
}

