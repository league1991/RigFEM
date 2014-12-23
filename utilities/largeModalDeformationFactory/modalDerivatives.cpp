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

#include "sparseMatrix.h"
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "StVKCubeABCD.h"
#include "StVKHessianTensor.h"
#include "StVKElementABCDLoader.h"
#include "insertRows.h"
#include "matrixIO.h"
#include "matrixPCA.h"
#include "computeStiffnessMatrixNullspace.h"
#include "largeModalDeformationFactory.h"
#include "sparseSolverAvailability.h"
#include "sparseSolvers.h"

// for faster computation, enable the -fopenmp -DUSE_OPENMP macro line in the Makefile-header file (see also documentation)

#ifdef USE_OPENMP
  #include <omp.h>
#endif

void MyFrame::OnLoadModalDerivatives(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Load modal derivatives"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.modalDeriv)|*.modalDeriv|All files(*.*)|*.*"), wxFD_OPEN /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString modalDerivativesFilename( dlg->GetPath());
    SaveCurrentWorkingDirectory(modalDerivativesFilename);
    if( !modalDerivativesFilename.empty() )
    {
      int m1,n1;
      double * newModalDerivatives = NULL;

      SetCursor(*wxHOURGLASS_CURSOR);
      const char * filename = modalDerivativesFilename.mb_str();
      int code = ReadMatrixFromDisk((char*)filename, 
        &m1, &n1, &newModalDerivatives);
      SetCursor(*wxSTANDARD_CURSOR);

      if (code != 0)
      {
        this->errMsg( _T("Loading error"),  
          _T("Unable to load modal derivatives from ") + modalDerivativesFilename );
        dlg->Destroy();
        return;
      }

      if (m1 != 3 * precomputationState.simulationMesh->getNumVertices())
      {
        this->errMsg( _T("Loading error"),  
          _T("The number of vertices in ") + modalDerivativesFilename + 
          _T(" does not match the simulation mesh."));
        free(newModalDerivatives);
        dlg->Destroy();
        return;
      }

      int targetNumDerivatives = (precomputationState.linearModalMatrix->Getr() - precomputationState.numRigidModes) * ( precomputationState.linearModalMatrix->Getr() - precomputationState.numRigidModes + 1) / 2;
      if (n1 != targetNumDerivatives)
      {
        char msg[4096];
        sprintf(msg, "The number of derivatives (%d; should be %d) is inconsistent with"
          " the number of linear modes (%d)" 
          " and the number of rigid modes (%d).",
          n1, targetNumDerivatives,
          precomputationState.linearModalMatrix->Getr(),
          precomputationState.numRigidModes);
        this->errMsg( _T("Loading error"), wxString(msg, wxConvUTF8));
        free(newModalDerivatives);
        dlg->Destroy();
        return;
      }

      // success
      delete(precomputationState.modalDerivativesMatrix);
      precomputationState.modalDerivativesMatrix = 
        new ModalMatrix(precomputationState.simulationMesh->getNumVertices(), 
        targetNumDerivatives, newModalDerivatives);
      free(newModalDerivatives);

      precomputationState.modalDerivativesAvailable = true;

      UpdateMenus();

      Refresh();
    }
  }

  dlg->Destroy();
}

void MyFrame::OnSaveModalDerivatives(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Save modal derivatives"), uiState.currentWorkingDirectory, _T(""), _T("Modal Matrix Files(*.modalDeriv)|*.modalDeriv|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);

  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString modalDerivativesFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(modalDerivativesFilename);
    if( !modalDerivativesFilename.empty() )
    {
      const char * filename = modalDerivativesFilename.mb_str();
      int code = WriteMatrixToDisk((char*)filename, 
        3 * precomputationState.modalDerivativesMatrix->Getn(), 
        precomputationState.modalDerivativesMatrix->Getr(), 
        precomputationState.modalDerivativesMatrix->GetMatrix());

      if (code != 0)
      {
        this->errMsg( _T("Saving error"),  
          _T("Unable to save modal derivatives to ") + modalDerivativesFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnExportModalDerivatives(wxCommandEvent& event)
{
  wxFileDialog *dlg = new wxFileDialog(this, _T("Export modal derivatives"), uiState.currentWorkingDirectory, _T(""), _T("Text Files(*.txt)|*.txt|All files(*.*)|*.*"), wxFD_SAVE /*| wxHIDE_READONLY*/, wxDefaultPosition);
  if ( dlg->ShowModal() == wxID_OK )
  {
    wxString modalDerivativesFilename( dlg->GetPath() );
    SaveCurrentWorkingDirectory(modalDerivativesFilename);
    if( !modalDerivativesFilename.empty() )
    {
      const char * filename = modalDerivativesFilename.mb_str();
      FILE * fout = fopen((char*)filename, "w");
      if (fout)
      {
        double * U = precomputationState.modalDerivativesMatrix->GetMatrix();
        int n = precomputationState.modalDerivativesMatrix->Getn();
        int r = precomputationState.modalDerivativesMatrix->Getr();
        fprintf(fout,"%d\n%d\n", 3*n, r);
        for(int i=0; i<3*n; i++)
        {
          for(int j=0; j<r; j++)
            fprintf(fout, "%.15f ", U[ELT(3*n,i,j)]);
          fprintf(fout,"\n");
        }
        fclose(fout);
      }
      else
      {
        this->errMsg( _T("Exporting error"),  
          _T("Unable to export modal derivatives to ") + modalDerivativesFilename );
        dlg->Destroy();
        return;
      }
    }
  }

  dlg->Destroy();
}

void MyFrame::OnComputeModalDerivatives(wxCommandEvent & event)
{
  if (!precomputationState.linearModesAvailable)
  {
    this->errMsg( _T("Cannot compute modal derivatives"),  
      _T("Linear modes are not available.\n"
      ) );
    return;
  }

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
    return;
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
}

void MyFrame::ComputeModalDerivatives(int * code, double ** modalDerivatives)
{
  *code = 0;

  // create stiffness matrix
  StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(precomputationState.simulationMesh);
  StVKInternalForces * internalForces = 
    new StVKInternalForces(precomputationState.simulationMesh, precomputedIntegrals); 

  // create stiffness matrix
  int n3 = 3 * precomputationState.simulationMesh->getNumVertices();
  SparseMatrix * stiffnessMatrix;
  StVKStiffnessMatrix * stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
  stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
  double * zero = (double*) calloc(n3, sizeof(double));
  stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);
  free(zero);

  // now, the stiffness matrix is computed
  // constrain the degrees of freedom
  int numConstrainedVertices = (int) (precomputationState.fixedVertices.size());
  int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
  set<int> :: iterator iter;
  int i = 0;
  for(iter = precomputationState.fixedVertices.begin(); 
      iter != precomputationState.fixedVertices.end(); 
      iter++)
  {
    constrainedDOFs[3*i+0] = 3 * (*iter) + 1;
    constrainedDOFs[3*i+1] = 3 * (*iter) + 2;
    constrainedDOFs[3*i+2] = 3 * (*iter) + 3;
    i++;
  }

  int oneIndexed = 1;
  stiffnessMatrix->RemoveRowsColumns(
    3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  int numRetainedDOFs = stiffnessMatrix->Getn();

  // generate rhs side
  bool computeHessianAtZero = (precomputationState.simulationMesh->getNumElements() < 5000);

  if (computeHessianAtZero)
    printf("Hessian at zero will be computed explicitly.\n");
  else
    printf("Hessian at zero will not be computed explicitly.\n");

  StVKHessianTensor * stVKStiffnessHessian = new 
    StVKHessianTensor(stiffnessMatrixClass);

  int numUsedLinearModes = precomputationState.rLin - precomputationState.numRigidModes;
  precomputationState.numDeriv = numUsedLinearModes * (numUsedLinearModes + 1) / 2;
  double * rhs = (double*) malloc (sizeof(double) * n3 * precomputationState.numDeriv);
  if (!rhs)
  {
    printf("Error: could not allocate space for all modal derivatives.\n");
    *code = 1;
    delete(precomputedIntegrals);
    delete(stiffnessMatrixClass);
    delete(internalForces);
    return;
  }

  if (computeHessianAtZero)
  {
    // compute hessian at zero
    if (stVKStiffnessHessian->ComputeHessianAtZero() != 0)
    {
      printf("Error: failed to evaluate the Hessian at the origin.\n");
      *code = 1;
      delete(precomputedIntegrals);
      delete(stiffnessMatrixClass);
      delete(internalForces);
      return;
    }
  }

  printf("Preparing to compute %d modal derivatives...\n", precomputationState.numDeriv);

  double * Ulin = precomputationState.linearModalMatrix->GetMatrix();
  if (computeHessianAtZero)
  {
    printf("Using the high-memory version.\n");
    int pos = 0;
    for(int i=0; i<numUsedLinearModes; i++)
    {
      printf("%d: ",i);fflush(NULL);
      for(int j=i; j<numUsedLinearModes; j++)
      {
        printf("%d ",j);fflush(NULL);
        stVKStiffnessHessian->EvaluateHessianQuadraticForm(
          &Ulin[n3*(precomputationState.numRigidModes + i)], &Ulin[n3*(precomputationState.numRigidModes + j)], &rhs[ELT(n3,0,pos)]);
        
        for(int k=0; k<n3; k++) //multiply by -1
          rhs[ELT(n3,k,pos)] *= -1.0;

        pos++;
      }
      printf("\n");
    }
  }
  else
  {
    printf("Using the low-memory version.\n");
    stVKStiffnessHessian->EvaluateHessianQuadraticFormDirectAll(
      Ulin,precomputationState.rLin,rhs,precomputationState.numRigidModes);
    
    if (n3*precomputationState.numDeriv < 0)
    {
      printf("Error: data too large to be indexed with the word size of your machine.\n");
      *code = 2;
      delete(stVKStiffnessHessian);
      delete(stiffnessMatrix);
      free(rhs);
      delete(precomputedIntegrals);
      delete(stiffnessMatrixClass);
      delete(internalForces);
      return;
    }

    /*
    if ((n3 > 200000) || (precomputationState.numDeriv > 1000))
    {
      printf("Warning: size of data %d might be too large to be indexed with the word size of your machine.\n",n3*precomputationState.numDeriv);
    }
    */

    //multiply by -1
    for(int i=0; i<n3*precomputationState.numDeriv; i++)
      rhs[i] *= -1.0;
  }

  printf("Right-hand sides for modal derivatives computed.\n");fflush(NULL);

  delete(stVKStiffnessHessian);
  delete(precomputedIntegrals);
  delete(stiffnessMatrixClass);
  delete(internalForces);

  // create mass matrix
  SparseMatrix * massMatrix;
  GenerateMassMatrix::computeMassMatrix(precomputationState.simulationMesh, &massMatrix, true);

  if (precomputationState.numRigidModes < 6)
  {
    double * buffer0 = (double*) malloc (sizeof(double) * n3);
    for(int i=0; i<precomputationState.numDeriv; i++)
    {
      massMatrix->MultiplyVector(&rhs[ELT(n3,0,i)], buffer0);
      for(int j=0; j<precomputationState.numRigidModes; j++)
      {
        // rhs -= <rhs, rigid mode j> * rigid mode j
        // rigid modes are mass-orthonormal
        double dotp = 0.0;
        for(int k=0; k<n3; k++)
          dotp += buffer0[k] * Ulin[ELT(n3,k,j)]; 
        for(int k=0; k<n3; k++)
          rhs[ELT(n3,k,i)] -= dotp * Ulin[ELT(n3,k,j)];
      }
    }
    free(buffer0);
  }
  else
  {
    RemoveSixRigidModes(precomputationState.numDeriv, rhs);
  }

  // constrain rhs
  double * rhsConstrained = (double*) malloc (sizeof(double) * numRetainedDOFs * precomputationState.numDeriv);
  for(int i=0; i<precomputationState.numDeriv; i++)
    RemoveRows(n3, &rhsConstrained[numRetainedDOFs*i], 
      &rhs[n3*i], 3 * numConstrainedVertices, constrainedDOFs, oneIndexed);

  free(rhs);

  // make room for (uninflated) derivatives
  double * modalDerivativesConstrained = (double*) malloc (sizeof(double) * precomputationState.numDeriv * numRetainedDOFs);

  // solve K * modesTemp = rhs
  printf("Factoring the %d x %d stiffness matrix...\n", numRetainedDOFs, numRetainedDOFs);
  //SPOOLESSolver * solver = new SPOOLESSolver(stiffnessMatrix);

  LinearSolver * solver;

  #ifdef PARDISO_SOLVER_IS_AVAILABLE
    int positiveDefinite = 0;
    int directIterative = 0;
    int numThreads = wxThread::GetCPUCount();
    PardisoSolver * pardisoSolver = new PardisoSolver(stiffnessMatrix, numThreads, positiveDefinite, directIterative);
    pardisoSolver->ComputeCholeskyDecomposition(stiffnessMatrix);
    solver = pardisoSolver;
  #elif defined(SPOOLES_SOLVER_IS_AVAILABLE)
    int numThreads = wxThread::GetCPUCount();
    if (numThreads > 1)
      solver = new SPOOLESSolverMT(stiffnessMatrix, numThreads);
    else
      solver = new SPOOLESSolver(stiffnessMatrix);
  #else
    solver = new CGSolver(stiffnessMatrix);
  #endif

  #ifdef USE_OPENMP
    #pragma omp parallel for
  #endif
  for(int i=0; i< precomputationState.numDeriv; i++)
  {
    printf("Solving for derivative #%d out of %d.\n", i + 1, precomputationState.numDeriv); fflush(NULL);
    solver->SolveLinearSystem(&modalDerivativesConstrained[ELT(numRetainedDOFs, 0, i)], &rhsConstrained[ELT(numRetainedDOFs, 0, i)]);
  }

  free(rhsConstrained);
  delete(solver);
  delete(stiffnessMatrix);

  *modalDerivatives = (double*) malloc (sizeof(double) * precomputationState.numDeriv * n3);

  // insert zero rows into the computed derivatives
  for(int i=0; i<precomputationState.numDeriv; i++)
  {
    InsertRows(n3, &modalDerivativesConstrained[numRetainedDOFs*i], 
      &((*modalDerivatives)[n3*i]), 
      3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
  }

  free(modalDerivativesConstrained);

  // remove rigid modes from modal derivatives
  if (precomputationState.numRigidModes > 0)
    printf("Removing rigid modes from modal derivatives...\n");fflush(NULL);

  if (precomputationState.numRigidModes < 6)
  {
    double * buffer1 = (double*) malloc (sizeof(double) * n3);
    for(int i=0; i<precomputationState.numDeriv; i++)
    {
      massMatrix->MultiplyVector(&((*modalDerivatives)[n3 * i]), buffer1);
      for(int j=0; j<precomputationState.numRigidModes; j++)
      {
        // rhs -= <rhs, rigid mode j> * rigid mode j
        // rigid modes are mass-orthonormal
        double dotp = 0.0;
        for(int k=0; k<n3; k++)
          dotp += buffer1[k] * Ulin[ELT(n3,k,j)]; 
        for(int k=0; k<n3; k++)
          (*modalDerivatives)[ELT(n3,k,i)] -= dotp * Ulin[ELT(n3,k,j)];
      }
    }
    free(buffer1);
  }
  else
  {
    RemoveSixRigidModes(precomputationState.numDeriv, *modalDerivatives);
  }

  //printf("Mass-normalizing modal derivatives...\n");fflush(NULL);

  // mass-normalize modal derivatives
  //for(int i=0; i < precomputationState.numDeriv; i++)
    //massMatrix->NormalizeVector(&((*modalDerivatives)[n3 * i]));

  delete(massMatrix);

  free(constrainedDOFs);
}

void MyFrame::RemoveSixRigidModes(int numVectors, double * x)
{
  int n3 = 3 * (precomputationState.simulationMesh)->getNumVertices();

  // remove six rigid modes from rhs
  double * nullspace6 = (double*) malloc (sizeof(double) * n3 * 6);
  double * defoPos6 = (double*) malloc (sizeof(double) * n3);
  for(int i=0; i<n3/3; i++)
  {
    Vec3d restPos = *((precomputationState.simulationMesh)->getVertex(i));
    for(int j=0; j<3; j++)
      defoPos6[3*i+j] = restPos[j];
  }

  ComputeStiffnessMatrixNullspace::ComputeNullspace(n3 / 3, defoPos6, nullspace6, 1, 1);
  free(defoPos6);

  for(int i=0; i<numVectors; i++)
  {
    ComputeStiffnessMatrixNullspace::RemoveNullspaceComponent(n3 / 3, 6, nullspace6, &x[ELT(n3,0,i)]);
  }

  free(nullspace6);
}

