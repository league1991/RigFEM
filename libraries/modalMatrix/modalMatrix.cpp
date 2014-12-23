/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "modalMatrix" library , Copyright (C) 2007 CMU                        *
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

#include <stdlib.h>
#include <string.h>
#include "lapack-headers.h"
#include "modalMatrix.h"

ModalMatrix::ModalMatrix(int n, int r, double * U, int flag)
{
  this->n = n;
  this->r = r;
  this->flag = flag;

  if (flag == 0)
  {
    this->U = (double*) malloc (sizeof(double) * 3 * n * r);
    memcpy(this->U,U,sizeof(double) * 3 * n * r);
  }
  else
  {
    this->U = U;
  }
}

ModalMatrix::~ModalMatrix()
{
  if (flag == 0)
    free(U);
}

void ModalMatrix::ProjectSingleVertex(int vertex, double vx, double vy, double vz, double * vreduced) 
{
  for (int j=0; j<r; j++) // over all columns of U
  {
    vreduced[j] = U[ELT(3*n,3*vertex+0,j)] * vx +
                  U[ELT(3*n,3*vertex+1,j)] * vy +
                  U[ELT(3*n,3*vertex+2,j)] * vz;
  }
}

void ModalMatrix::AddProjectSingleVertex(int vertex, double vx, double vy, double vz, double * vreduced) 
{
  for (int j=0; j<r; j++) // over all columns of U
  {
    vreduced[j] += 
      U[ELT(3*n,3*vertex+0,j)] * vx +
      U[ELT(3*n,3*vertex+1,j)] * vy +
      U[ELT(3*n,3*vertex+2,j)] * vz;
  }
}

void ModalMatrix::ProjectVector(double * v, double * vreduced) 
{
  // has to make inner product of vector f will all the columns of U
  // i.e. multiply U^T * f = q
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasTrans;
  int M = 3*n;
  int N = r;
  double alpha = 1;
  double * a = U;
  int lda = 3*n;
  double * x = v;
  int incx = 1;
  double beta = 0;
  double * y = vreduced; 
  int incy = 1;

  cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

void ModalMatrix::AddProjectVector(double * v, double * vreduced) 
{
  // has to make inner product of vector f will all the columns of U
  // i.e. multiply U^T * f = q
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasTrans;
  int M = 3*n;
  int N = r;
  double alpha = 1;
  double * a = U;
  int lda = 3*n;
  double * x = v;
  int incx = 1;
  double beta = 1;
  double * y = vreduced; 
  int incy = 1;

  cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

void ModalMatrix::ProjectSparseVector(int numSparseEntries, double * sparseVector, int * sparseVectorIndices, double * vreduced)
{
  // has to make inner product of vector f will all the columns of U
  int i,j;
  for (j=0; j<r; j++) // over all columns of U
  {
    // dot product of column j of U with vector f
    vreduced[j] = 0;
    for (i=0; i<numSparseEntries; i++)
      vreduced[j] += U[ELT(3*n,sparseVectorIndices[i],j)] * sparseVector[i];
  }
}

void ModalMatrix::AddProjectSparseVector(int numSparseEntries, double * sparseVector, int * sparseVectorIndices, double * vreduced)
{
  // has to make inner product of vector f will all the columns of U
  int i,j;
  for (j=0; j<r; j++) // over all columns of U
  {
    // dot product of column j of U with vector f
    for (i=0; i<numSparseEntries; i++)
      vreduced[j] += U[ELT(3*n,sparseVectorIndices[i],j)] * sparseVector[i];
  }
}

void ModalMatrix::ProjectMatrix(int numColumns, double * matrix, double * matrixReduced)
{
  // multiply U^T * matrix = matrixReduced

  CBLAS_ORDER order= CblasColMajor;
  int M = r;
  int N = numColumns;
  int K = 3 * n;
  double alpha = 1.0;
  double * a = U;
  int lda = K;
  double * b = matrix;
  int ldb = 3 * n;
  double beta = 0.0;
  double * c = matrixReduced; 
  int ldc = r;

  cblas_dgemm(order, CblasTrans, CblasNoTrans,
      M, N, K, alpha, a, lda, b, ldb,
      beta, c, ldc);
}

void ModalMatrix::AssembleMatrix(int numColumns, double * qMatrix, double * uMatrix)
{
  for(int i=0; i<numColumns; i++)
    AssembleVector(&qMatrix[ELT(r,0,i)], &uMatrix[ELT(3*n,0,i)]);
}

void ModalMatrix::AssembleVector(double * q, double * u) //u = U * q;
{
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasNoTrans;
  int M = 3*n;
  int N = r;
  double alpha = 1;
  double * a = U;
  int lda = 3*n;
  double * x = q;
  int incx = 1;
  double beta = 0;
  double * y = u; 
  int incy = 1;

  cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

void ModalMatrix::AddAssembleVector(double * q, double * u) //u = U * q;
{
  CBLAS_ORDER     order= CblasColMajor;
  CBLAS_TRANSPOSE trans= CblasNoTrans;
  int M = 3*n;
  int N = r;
  double alpha = 1;
  double * a = U;
  int lda = 3*n;
  double * x = q;
  int incx = 1;
  double beta = 1.0;
  double * y = u; 
  int incy = 1;

  cblas_dgemv(order, trans, M, N, alpha, a, lda, x, incx, beta, y, incy);
}

