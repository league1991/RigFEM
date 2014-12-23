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

#ifndef _MODALMATRIX_H_
#define _MODALMATRIX_H_

#include "matrixMacros.h"

#define MODAL_MATRIX_INTERNAL_COPY 0
#define MODAL_MATRIX_NO_INTERNAL_COPY 1

class ModalMatrix
{
public:
  // n is num vertices, matrix U is 3n x r
  ModalMatrix(int n, int r, double * U, int flag=0); // flag = 0: makes an internal copy of U; flag != 0: does not make an internal copy of U
  ~ModalMatrix();

  inline int Getr() { return r; }
  inline int Getn() { return n; }

  inline double * GetMatrix() { return U; }

  // computes: vreduced = U^T * v
  void ProjectVector(double * v, double * vreduced);
  void ProjectSparseVector(int numSparseEntries, double * sparseVector, int * sparseVectorIndices, double * vreduced);
  void ProjectSingleVertex(int vertex, double vx, double vy, double vz, double * vreduced);

  // computes matrixReduced = U^T * matrix; matrix is 3n x numColumns
  void ProjectMatrix(int numColumns, double * matrix, double * matrixReduced);

  // vreduced += U^T * v
  void AddProjectVector(double * v, double * vreduced);
  void AddProjectSparseVector(int numSparseEntries, double * sparseVector, int * sparseVectorIndices, double * vreduced);
  void AddProjectSingleVertex(int vertex, double vx, double vy, double vz, double * vreduced);

  // computes u = U * q
  void AssembleVector(double * q, double * u);
  void AssembleMatrix(int numColumns, double * qMatrix, double * uMatrix);
  inline void AssembleSingleVertex(int vertex, double * q, double * ux, double * uy, double * uz);

  // u += U * q
  void AddAssembleVector(double * q, double * u);
  inline void AddAssembleSingleVertex(int vertex, double * q, double * ux, double * uy, double * uz);

protected:

  double * U; // pointer to the deformation basis

  int r; // number of columns
  int n; // number of vertices

  int flag;
};

// constructs the deformation of vertex 'vertex', given q
// result goes into (ux,uy,uz)
inline void ModalMatrix::AssembleSingleVertex(int vertex, double * q, double * ux, double * uy, double * uz)
{
  double regx = 0;
  double regy = 0;
  double regz = 0;

  double * pos = &U[ELT(3*n,3*vertex,0)];
  int n3 = 3*n;

  for (int j=0; j<r; j++) // over all columns of U
  {
    regx += pos[0] * q[j];
    regy += pos[1] * q[j];
    regz += pos[2] * q[j];
    pos += n3;
  }

  *ux = regx;
  *uy = regy;
  *uz = regz;
}

inline void ModalMatrix::AddAssembleSingleVertex(int vertex, double * q, double * ux, double * uy, double * uz)
{
  double regx = 0;
  double regy = 0;
  double regz = 0;

  double * pos = &U[ELT(3*n,3*vertex,0)];
  int n3 = 3*n;

  for (int j=0; j<r; j++) // over all columns of U
  {
    regx += pos[0] * q[j];
    regy += pos[1] * q[j];
    regz += pos[2] * q[j];
    pos += n3;
  }

  *ux += regx;
  *uy += regy;
  *uz += regz;
}

#endif

