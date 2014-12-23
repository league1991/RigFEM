/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "elasticForceModel" library , Copyright (C) 2007 CMU, 2009 MIT,       *
 *                                                       2014 USC        *
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
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "matrixMacros.h"
#include "performanceCounter.h"
#include "reducedMassSpringSystemForceModel.h"

#if USE_MKL_SPARSE_BLAS
  #include "mkl_spblas.h"
  #include "mkl.h"
#endif

ReducedMassSpringSystemForceModel::ReducedMassSpringSystemForceModel(
                         MassSpringSystem * massSpringSystem_,
                         ModalMatrix * modalMatrix_): ReducedForceModel(), massSpringSystem(massSpringSystem_), modalMatrix(modalMatrix_)
{
  massSpringSystem->GetStiffnessMatrixTopology(&sparseMatrix);

  n = modalMatrix->Getn();
  r = modalMatrix->Getr();
  U = modalMatrix->GetMatrix();

  // allocate buffers
  u = (double*) malloc (sizeof(double) * 3 * n);
  bufferVector = (double*) malloc (sizeof(double) * 3 * n);
  bufferMatrix = (double*) malloc (sizeof(double) * 3 * n * r);

  #if USE_MKL_SPARSE_BLAS
    int nnz = sparseMatrix->GetNumUpperTriangleEntries();
    csr_values = (double*) malloc (sizeof(double) * nnz);
    csr_columns = (int*) malloc (sizeof(int) * nnz);
    csr_pointerB = (int*) malloc (sizeof(int) * sparseMatrix->GetNumRows());
    csr_pointerE = (int*) malloc (sizeof(int) * sparseMatrix->GetNumRows());
  #endif
}

ReducedMassSpringSystemForceModel::~ReducedMassSpringSystemForceModel()
{
  free(bufferVector);
  free(bufferMatrix);
  free(u);
  delete(sparseMatrix);
  #if USE_MKL_SPARSE_BLAS
    free(csr_values);
    free(csr_columns);
    free(csr_pointerB);
    free(csr_pointerE);
  #endif
}

void ReducedMassSpringSystemForceModel::GetInternalForce(double * q, double * internalForces)
{
  // construct u = U * q
  modalMatrix->AssembleVector(q, u);

  // evaluate forces
  massSpringSystem->ComputeForce(u, bufferVector);

  // project forces
  modalMatrix->ProjectVector(bufferVector, internalForces);

  for(int i=0; i<r; i++)
    internalForces[i] *= -1;
}

void ReducedMassSpringSystemForceModel::GetTangentStiffnessMatrix(
  double * q, double * tangentStiffnessMatrix)
{
  // construct u = U * q
  modalMatrix->AssembleVector(q, u);
  GetTangentStiffnessMatrixHelper(tangentStiffnessMatrix);
}

void ReducedMassSpringSystemForceModel::GetForceAndMatrix(double * q, 
    double * internalForces, double * tangentStiffnessMatrix)
{
  GetInternalForce(q, internalForces);
  GetTangentStiffnessMatrixHelper(tangentStiffnessMatrix);
}

void ReducedMassSpringSystemForceModel::GetTangentStiffnessMatrixHelper(
  double * tangentStiffnessMatrix)
{
  // evaluate stiffness matrix
  //PerformanceCounter counter;
  massSpringSystem->ComputeStiffnessMatrix(u, sparseMatrix);
  //counter.StopCounter();
  //printf("counter: %G\n", counter.GetElapsedTime());

  // project matrix
  #if USE_MKL_SPARSE_BLAS
    mkl_set_num_threads(8);
    //PerformanceCounter counter;

    int upperTriangleOnly=1;
    int oneIndexed=1;
    sparseMatrix->GenerateCompressedRowMajorFormat_four_array(csr_values, csr_columns, csr_pointerB, csr_pointerE, upperTriangleOnly, oneIndexed); 

    char transa = 'N';
    int m = sparseMatrix->GetNumRows();
    int n = r;
    int k = m;
    double alpha = 1.0;
    char matdescra[7] = "SUNFXX";
    double * val = csr_values;
    int * indx = csr_columns;
    int * pntrb = csr_pointerB;
    int * pntre = csr_pointerE;
    double * b = U;
    int ldb = m;
    double beta = 0.0;
    double * c = bufferMatrix;
    int ldc = m;
    mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, val, indx, pntrb, pntre,
      b, &ldb, &beta, c, &ldc);

    //counter.StopCounter();
    //printf("counter: %G\n", counter.GetElapsedTime());
  #else
    for(int i=0; i<r; i++)
      sparseMatrix->MultiplyVector(&U[ELT(3*n,0,i)], &bufferMatrix[ELT(3*n,0,i)]);
  #endif

  modalMatrix->ProjectMatrix(r, bufferMatrix, tangentStiffnessMatrix);

  int r2 = r*r;
  for(int i=0; i<r2; i++)
    tangentStiffnessMatrix[i] *= -1;
}

