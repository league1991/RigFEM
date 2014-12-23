/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "corotational linear FEM" library , Copyright (C) 2014 USC            *
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

#ifndef _COROTATIONALLINEARFEM_H_
#define _COROTATIONALLINEARFEM_H_

/*
  Corotational linear FEM deformable model.

  This class implements the deformable model described in the following paper:

  M. Mueller, M. Gross: Interactive Virtual Materials.
  In Proc. of Graphics Interface 2004 (2004), pp. 239–246.

  In [Mueller 2004], the tangent stiffness matrix is approximate (warp=1). 
  This class can also compute the exact tangent stiffness matrix (warp=2).
  The implementation is described in:
  J. Barbic: Exact Corotational Linear FEM Stiffness Matrix, Technical Report, USC, 2012

  It is also possible to turn warping off (warp=0). This gives fast linear FEM dynamics,
  but large deformations are not well-represented.
*/

#include "tetMesh.h"
#include "sparseMatrix.h"

class CorotationalLinearFEM
{
public:

  // initializes corotational linear FEM
  // input: tetMesh
  CorotationalLinearFEM(TetMesh * tetMesh);
  virtual ~CorotationalLinearFEM();

  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); // returns a zero matrix containing the locations of non-zero elements in the stiffness matrix

  // computes the internal forces and (warped) stiffness matrix for the entire mesh
  // vertex displacements (input) and internal forces (output) must be (pre-allocated) vectors of length 3 * numVertices
  // the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
  // i.e., the computed internal forces are *negatives* of the actual physical internal forces acting on the material
  // warp:
  //   0: no warping (linear FEM)
  //   1: stiffness warping (corotational linear FEM with approximate stiffness matrix) [Mueller 2004]
  //   2: corotational linear FEM with exact tangent stiffness matrix (see the technical report [Barbic 2012])
  virtual void ComputeForceAndStiffnessMatrix(double * vertexDisplacements, double * internalForces, SparseMatrix * stiffnessMatrix, int warp=1);

  // this routine is same as above, except that it only traverses elements from elementLo <= element <= elementHi - 1
  void ComputeForceAndStiffnessMatrixOfSubmesh(double * vertexDisplacements, double * internalForces, SparseMatrix * stiffnessMatrix, int warp, int elementLo, int elementHi);

  inline TetMesh * GetTetMesh() { return tetMesh; }

protected:
  int numVertices;
  TetMesh * tetMesh;
  double * undeformedPositions;
  double ** MInverse;
  double ** KElementUndeformed;

  void WarpMatrix(double * K, double * R, double * RK, double * RKRT);
  void inverse3x3(double * A, double * AInv); // inverse of a row-major 3x3 matrix
  void inverse4x4(double * A, double * AInv); // inverse of a row-major 4x4 matrix

  // acceleration indices
  int ** rowIndices;
  int ** columnIndices;
  void ClearRowColumnIndices();
  void BuildRowColumnIndices(SparseMatrix * sparseMatrix);
};

#endif

