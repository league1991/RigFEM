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

#include "reducedMassSpringSystemForceModelWithHessian.h"

#if USE_MKL_SPARSE_BLAS
  #include "mkl.h"
  #include "mkl_spblas.h"
#endif

ReducedMassSpringSystemForceModelWithHessian::ReducedMassSpringSystemForceModelWithHessian(MassSpringSystem * massSpringSystem, ModalMatrix * modalMatrix): ReducedMassSpringSystemForceModel(massSpringSystem, modalMatrix), ReducedForceModelWithHessian()
{
  bufferTangentMatrix = (double*) malloc (sizeof(double) * r * r);
}

ReducedMassSpringSystemForceModelWithHessian::~ReducedMassSpringSystemForceModelWithHessian()
{
  free(bufferTangentMatrix);
}

void ReducedMassSpringSystemForceModelWithHessian::GetTangentHessianTensor(double * q, double * tangentHessianTensor)
{
  modalMatrix->AssembleVector(q,u);
  for(int i=0; i<r; i++)
  {
    massSpringSystem->ComputeStiffnessMatrixCorrection(u, &U[ELT(3*n,0,i)], sparseMatrix);
    // project matrix
    for(int j=0; j<r; j++)
      sparseMatrix->MultiplyVector(&U[ELT(3*n,0,j)], &bufferMatrix[ELT(3*n,0,j)]);
    modalMatrix->ProjectMatrix(r, bufferMatrix, bufferTangentMatrix);

    int r2 = r*r;
    for(int j=0; j<r2; j++)
      bufferTangentMatrix[j] *= -1;

    // write it in place
    for(int l=r-1; l>=0; l--)
      for(int k=r-1; k>=l; k--)
      {
        int lowerTrianglePos = l * r - (l-1) * l / 2 + (k-l);
        tangentHessianTensor[ELT(r,i,lowerTrianglePos)] = bufferTangentMatrix[ELT(r,k,l)];
      }
  }
}

void ReducedMassSpringSystemForceModelWithHessian::GetStiffnessMatrixCorrection(double * q, double * dq, double * dK)
{
  modalMatrix->AssembleVector(q, u);
  modalMatrix->AssembleVector(dq, bufferVector);
  massSpringSystem->ComputeStiffnessMatrixCorrection(u, bufferVector, sparseMatrix);

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
    for(int j=0; j<r; j++)
      sparseMatrix->MultiplyVector(&U[ELT(3*n,0,j)], &bufferMatrix[ELT(3*n,0,j)]);
  #endif

  modalMatrix->ProjectMatrix(r, bufferMatrix, dK);

  int r2 = r*r;
  for(int j=0; j<r2; j++)
    dK[j] *= -1;
}

