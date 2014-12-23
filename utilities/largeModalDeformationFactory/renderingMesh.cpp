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

// manage the rendering mesh

#include "generateSurfaceMesh.h"
#include "largeModalDeformationFactory.h"

void MyFrame::OnLoadMesh(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load mesh"), uiState.currentWorkingDirectory, _T(""), _T("Obj Files(*.obj)|*.obj|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition); 

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString objMeshname(dlg->GetPath());
    if( !objMeshname.empty() )
    {
      ObjMesh * newObjMesh;
      SetCursor(*wxHOURGLASS_CURSOR);
      try
      {
        newObjMesh = new ObjMesh(string(objMeshname.mb_str()));
      }
      catch( ObjMeshException myException )
      {
        SetCursor(*wxSTANDARD_CURSOR);
        this->errMsg( _T("Loading error"),  wxString(myException.getReason().c_str(), wxConvUTF8) );
        dlg->Destroy();
        return;
      }

      SetCursor(*wxSTANDARD_CURSOR);

      precomputationState.renderingMeshAvailable = true;
      DeleteInterpolationData();
      UpdateMenus();

      delete(precomputationState.renderingMesh); // deallocate any old objMesh
      precomputationState.renderingMesh = newObjMesh;
      myGLCanvas->UpdateRenderingMeshRenderData();
      SelectView(UIState::VIEW_RENDERING_MESH);

      uiState.projectName = objMeshname.AfterLast('/').AfterLast('\\').BeforeLast('.');
      uiState.renderingMeshFilename = objMeshname;
    }
  }
  dlg->Destroy();
}

void MyFrame::CreateRenderingMeshFromSimulationMesh()
{
  delete(precomputationState.renderingMesh); // deallocate any old objMesh
  precomputationState.renderingMesh = GenerateSurfaceMesh::ComputeMesh(precomputationState.simulationMesh);
  precomputationState.renderingMeshAvailable = true;

  DeleteInterpolationData();
  UpdateMenus();

  int updateCamera = 0;
  myGLCanvas->UpdateRenderingMeshRenderData(updateCamera);
}

void MyFrame::DeleteRenderingMesh()
{
  delete(precomputationState.renderingMesh); // deallocate any old objMesh
  precomputationState.renderingMesh = NULL;
  precomputationState.renderingMeshAvailable = false;
  DeleteInterpolationData();
  UpdateMenus();
  myGLCanvas->DeleteRenderingMeshRenderData();
}

void MyFrame::OnRemoveTriangleMesh(wxCommandEvent& event)
{
  DeleteRenderingMesh();
}

