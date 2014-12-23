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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "corotationalLinearFEM.h"
#include "polarDecomposition.h"
#include "matrixMultiplyMacros.h"
#include "mat3d.h"
#include "volumetricMeshENuMaterial.h"
#include "volumetricMeshOrthotropicMaterial.h"

CorotationalLinearFEM::CorotationalLinearFEM(TetMesh * tetMesh_) : tetMesh(tetMesh_) 
{
  numVertices = tetMesh->getNumVertices();

  // store the undeformed positions
  undeformedPositions = (double*) malloc (sizeof(double) * 3 * numVertices);
  for(int i=0; i < numVertices; i++)
  {
    Vec3d * v = tetMesh->getVertex(i);
    for(int j=0; j<3; j++)
      undeformedPositions[3*i+j] = (*v)[j];
  }

  int numElements = tetMesh->getNumElements();

  MInverse = (double**) malloc (sizeof(double*) * numElements);
  for(int el = 0; el < numElements; el++)
  {
    // get the integer indices of the tet vertices
    int vtxIndex[4];
    for(int vtx=0; vtx<4; vtx++)
      vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);
    /*
       Form matrix: 
       M = [ v0   v1   v2   v3 ]
           [  1    1    1    1 ]
    */
    double M[16]; // row-major
    for(int vtx=0; vtx<4; vtx++)
      for(int dim=0; dim<3; dim++)
        M[4 * dim + vtx] = undeformedPositions[3 * vtxIndex[vtx] + dim];
    M[12] = M[13] = M[14] = M[15] = 1.0;

    // invert M and cache inverse (see [Mueller 2004])
    MInverse[el] = (double*) malloc (sizeof(double) * 16);
    inverse4x4(M, MInverse[el]);
  }

  // build acceleration indices for fast writing to the global stiffness matrix
  SparseMatrix * sparseMatrix;
  GetStiffnessMatrixTopology(&sparseMatrix);
  BuildRowColumnIndices(sparseMatrix);
  delete(sparseMatrix);

  // compute stiffness matrices for all the elements in the undeformed configuration
  KElementUndeformed = (double**) malloc (sizeof(double*) * numElements);
  for (int el = 0; el < numElements; el++)
  {
    double * MInv = MInverse[el];

    // Form stiffness matrix of the element in the undeformed configuration.
    // The procedure below is standard in FEM solid mechanics.
    // This code implements the equations given in Ahmed A. Shabana: Theory of Vibration, Volume II: Discrete and Continuous Systems, Springer--Verlag, New York, NY, 1990.

    double B[72] = 
      { MInv[0], 0, 0, MInv[4], 0, 0, MInv[8], 0, 0, MInv[12], 0, 0,
        0, MInv[1], 0, 0, MInv[5], 0, 0, MInv[9], 0, 0, MInv[13], 0,
	0, 0, MInv[2], 0, 0, MInv[6], 0, 0, MInv[10], 0, 0, MInv[14],
        MInv[1], MInv[0], 0, MInv[5], MInv[4], 0, MInv[9], MInv[8], 
        0, MInv[13], MInv[12], 0, 0, MInv[2], MInv[1], 0, MInv[6], MInv[5], 
        0, MInv[10], MInv[9], 0, MInv[14], MInv[13], MInv[2], 0, MInv[0], 
        MInv[6], 0, MInv[4], MInv[10], 0, MInv[8], MInv[14], 0, MInv[12] };

    // compute elasticity stiffness tensor
    double E[36]; // 6 x 6 matrix, stored row-major (symmetric, so row vs column-major storage makes no difference anyway)
    VolumetricMesh::Material * material = tetMesh->getElementMaterial(el);

    // check if material is ENuMaterial (i.e., isotropic)
    VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
    if (eNuMaterial != NULL)
    {
      // material is isotropic, specified by E, nu
      // compute Lame coefficients
      double lambda = eNuMaterial->getLambda();
      double mu = eNuMaterial->getMu();

      double Et[36] = { lambda + 2 * mu, lambda, lambda, 0, 0, 0,
                        lambda, lambda + 2 * mu, lambda, 0, 0, 0,
                        lambda, lambda, lambda + 2 * mu, 0, 0, 0,
                        0, 0, 0, mu, 0, 0,
                        0, 0, 0, 0, mu, 0,
                        0, 0, 0, 0, 0, mu };

      memcpy(E, Et, sizeof(double) * 36);
    }
    else
    {
      // orthotropic material
      // we follow the following references:
      // Yijing Li and Jernej Barbic: Stable Orthotropic Materials, Symposium on Computer Animation 2014
      // http://en.wikipedia.org/wiki/Orthotropic_material
      // http://www.solidmechanics.org/text/Chapter3_2/Chapter3_2.htm

      // test if material is OrthotropicMaterial (i.e., orthotropic)
      VolumetricMesh::OrthotropicMaterial * orthotropicMaterial = downcastOrthotropicMaterial(material);
      if (orthotropicMaterial != NULL)
      {
        double E1 = orthotropicMaterial->getE1();
        double E2 = orthotropicMaterial->getE2();
        double E3 = orthotropicMaterial->getE3();
        double nu12 = orthotropicMaterial->getNu12();
        double nu23 = orthotropicMaterial->getNu23();
        double nu31 = orthotropicMaterial->getNu31();
        double G12 = orthotropicMaterial->getG12();
        double G23 = orthotropicMaterial->getG23();
        double G31 = orthotropicMaterial->getG31();

        double nu21 = nu12 * E2 / E1;
        double nu32 = nu23 * E3 / E2;
        double nu13 = nu31 * E1 / E3;
        
        double Y = 1.0 / (1.0 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 - 2.0 * nu21 * nu32 * nu13);

        double ELocal[36] = { E1 * (1.0 - nu23 * nu32) * Y, E1 * (nu21 + nu31 * nu23) * Y, E1 * (nu31 + nu21 * nu32) * Y, 0.0, 0.0, 0.0,
                              E1 * (nu21 + nu31 * nu23) * Y, E2 * (1.0 - nu13 * nu31) * Y, E2 * (nu32 + nu12 * nu31) * Y, 0.0, 0.0, 0.0,
                              E1 * (nu31 + nu21 * nu32) * Y, E2 * (nu32 + nu12 * nu31) * Y, E3 * (1.0 - nu12 * nu21) * Y, 0.0, 0.0, 0.0,
                              0, 0, 0, G12, 0, 0,
                              0, 0, 0, 0, G23, 0,
                              0, 0, 0, 0, 0, G31 };

        //memcpy(E, ELocal, sizeof(double) * 36); // debug

        double R[9]; // row-major
        orthotropicMaterial->getR(R);

        // rotate Elocal into the basis given by the columns of R
        #define Relt(i,j) (R[3*(i)+(j)])
        double rotator[36];
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
          {
            rotator[6 * i + j] = Relt(i,j) * Relt(i,j);
            rotator[6 * i + 3 + j] = 2.0 * Relt(i, j) * Relt(i, (j+1) % 3);
            rotator[6 * (i + 3) + j] = Relt(i, j) * Relt((i+1) % 3, j);
            rotator[6 * (i + 3) + 3 + j] = Relt(i, j) * Relt((i+1) % 3, (j+1) % 3) + Relt(i, (j+1) % 3) * Relt((i+1) % 3, j);
          }
        #undef Relt

        // debug
        //memset(rotator, 0, sizeof(double) * 36);
        //for(int i=0; i<6; i++)
          //rotator[6*i+i] = 1.0;

        // E = rotator * ELocal * rotator^T
        double buffer[36]; 
        memset(buffer, 0, sizeof(double) * 36);
        // buffer = ELocal * rotator^T
        for(int i=0; i<6; i++)
          for(int j=0; j<6; j++)
            for(int k=0; k<6; k++)
              buffer[6 * i + j] += ELocal[6 * i + k] * rotator[6 * j + k];

        // E = rotator * buffer
        memset(E, 0, sizeof(double) * 36);
        for(int i=0; i<6; i++)
          for(int j=0; j<6; j++)
            for(int k=0; k<6; k++)
              E[6 * i + j] += rotator[6 * i + k] * buffer[6 * k + j];

      }
      else
      {
        printf("Error: CorotationalLinearFEM: unknown material encounted in the mesh.\n");
        throw 1;
      }
    }

    // EB = E * B
    double EB[72];
    memset(EB, 0, sizeof(double) * 72);
    for (int i=0; i<6; i++)
      for (int j=0; j<12; j++)
	for (int k=0; k<6; k++)
	  EB[12 * i + j] += E[6 * i + k] * B[12 * k + j];
 
    // KElementUndeformed[el] = B^T * EB
    KElementUndeformed[el] = (double*) calloc (144, sizeof(double)); // element stiffness matrix
    for (int i=0; i<12; i++)
      for (int j=0; j<12; j++)
	for (int k=0; k<6; k++)
          KElementUndeformed[el][12 * i + j] += B[12 * k + i] * EB[12 * k + j];

    // KElementUndeformed[el] *= volume
    double volume = TetMesh::getTetVolume(tetMesh->getVertex(el,0), tetMesh->getVertex(el,1), tetMesh->getVertex(el,2), tetMesh->getVertex(el,3));
    for(int i=0; i<144; i++)
      KElementUndeformed[el][i] *= volume;
  }
}

CorotationalLinearFEM::~CorotationalLinearFEM()
{
  free(undeformedPositions);
  for(int el=0; el < tetMesh->getNumElements(); el++)
  {
    free(KElementUndeformed[el]);
    free(MInverse[el]);
  }
  free(KElementUndeformed);
  free(MInverse);

  ClearRowColumnIndices();
}

void CorotationalLinearFEM::GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology)
{
  SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);

  int numElements = tetMesh->getNumElements();
  for (int el=0; el < numElements; el++)
  {
    int vtxIndex[4];
    for(int vtx=0; vtx<4; vtx++)
      vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

    for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
      {
        // add 3x3 block corresponding to pair of vertices (i,j)
        for(int k=0; k<3; k++)
          for(int l=0; l<3; l++)
            emptyMatrix->AddEntry(3 * vtxIndex[i] + k, 3 * vtxIndex[j] + l, 0.0);
      }
  }

  *stiffnessMatrixTopology = new SparseMatrix(emptyMatrix);
  delete(emptyMatrix);
}

// compute RK = R * K and RKRT = R * K * R^T (block-wise)
// input: K, R
// output: RK, RKRT
void CorotationalLinearFEM::WarpMatrix(double * K, double * R, double * RK, double * RKRT)
{
  memset(RK, 0, sizeof(double) * 144);
  memset(RKRT, 0, sizeof(double) * 144);
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
    {
      // RK = R * K
      for(int k=0; k<3; k++)
         for(int l=0; l<3; l++)
           for(int m=0; m<3; m++)
             RK[12 * (3 * i + k) + (3 * j + l)] += R[3 * k + m] * K[12 * (3 * i + m) + (3 * j + l)];

      // RKRT = RK * R^T
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          for(int m=0; m<3; m++)
            RKRT[12 * (3 * i + k) + (3 * j + l)] += RK[12 * (3 * i + k) + (3 * j + m)] * R[3 * l + m];
    }
}

void CorotationalLinearFEM::ComputeForceAndStiffnessMatrix(double * u, double * f, SparseMatrix * stiffnessMatrix, int warp)
{
  ComputeForceAndStiffnessMatrixOfSubmesh(u, f, stiffnessMatrix, warp, 0, tetMesh->getNumElements());
}

void CorotationalLinearFEM::ComputeForceAndStiffnessMatrixOfSubmesh(double * u, double * f, SparseMatrix * stiffnessMatrix, int warp, int elementLo, int elementHi)
{
  // clear f to zero
  if (f != NULL)
    memset(f, 0, sizeof(double) * 3 * numVertices);

  // clear stiffness matrix to zero
  if (stiffnessMatrix != NULL)
    stiffnessMatrix->ResetToZero();

  for (int el=elementLo; el < elementHi; el++)
  {
    int vtxIndex[4];
    for (int vtx=0; vtx<4; vtx++)
      vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

    double KElement[144]; // element stiffness matrix, to be computed below; row-major

    if (warp > 0)
    {
      double P[16]; // the current world-coordinate positions (row-major)
      /*
         P = [ v0   v1   v2   v3 ]
             [  1    1    1    1 ]
      */
      // rows 1,2,3
      for(int i=0; i<3; i++)
        for(int j=0; j<4; j++)
          P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
      // row 4
      for(int j=0; j<4; j++)
        P[12 + j] = 1;

      // F = P * Inverse(M)
      double F[9]; // upper-left 3x3 block
      for(int i=0; i<3; i++) 
        for(int j=0; j<3; j++) 
        {
          F[3 * i + j] = 0;
          for(int k=0; k<4; k++)
            F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
	}

      double R[9]; // rotation (row-major)
      double S[9]; // symmetric (row-major)
      double tolerance = 1E-6;
      int forceRotation = 1;
      PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

      // RK = R * K
      // KElement = R * K * R^T
      double RK[144]; // row-major
      WarpMatrix(KElementUndeformed[el], R, RK, KElement);

      // f = RK (RT x - x0)
      double fElement[12];
      for(int i=0; i<12; i++)
      {
        fElement[i] = 0;
        for(int j=0; j<4; j++)
          for(int l=0; l<3; l++)
            fElement[i] += KElement[12 * i + 3 * j + l] * P[4 * l + j] - RK[12 * i + 3 * j + l] * undeformedPositions[3 * vtxIndex[j] + l];
      }

      // add fElement into the global f
      if (f != NULL)
      {
        for(int j=0; j<4; j++)
          for(int l=0; l<3; l++)
            f[3 * vtxIndex[j] + l] += fElement[3 * j + l];
      }

      // compute exact stiffness matrix
      if (warp == 2)
      {
        // compute G = (tr(S) I - S) R^T
        double G[9]; 
        double tr = S[0] + S[4] + S[8];
        double temp[9];
        for(int i=0; i<9; i++)
          temp[i] = -S[i];
        temp[0] += tr;
        temp[4] += tr;
        temp[8] += tr;
        // G = temp * R^T
        MATRIX_MULTIPLY3X3ABT(temp, R, G);

        double invG[9]; // invG = G^{-1}
        inverse3x3(G, invG);

        double rhs[27]; // 3 x 9 matrix (column-major)
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
          {
            double temp[9];
            for(int k=0; k<9; k++)
              temp[k] = 0.0;
            // copy i-th row of R into column j of temp      
            for(int k=0; k<3; k++)
              temp[3 * k + j] = R[3 * i + k];
            // extract the skew-symmetric part
            SKEW_PART(temp, &rhs[3 * (3 * i + j)]);
          }
        // must undo division by 2 from inside the SKEW_PART macro
        for(int i=0; i<27; i++)
          rhs[i] *= 2.0;

        // solve G * omega = rhs
        double omega[27]; // column-major
        for(int i=0; i<9; i++)
        {
          MATRIX_VECTOR_MULTIPLY3X3(invG, &rhs[3 * i], &omega[3 * i]);
        }

        double dRdF[81]; // each column is skew(omega) * R ; column-major
        for(int i=0; i<9; i++)
        {
          double skew[9];
          SKEW_MATRIX(&omega[3 * i], skew);
          MATRIX_MULTIPLY3X3(skew, R, &dRdF[9 * i]);
        }

        double B[3][3][9];
        // re-arrange dRdF into B, for easier dRdF * dFdx multiplication (to exploit sparsity of dFdx)
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
            for(int k=0; k<3; k++)
              for(int l=0; l<3; l++)
              {
                int row = 3 * i + k;
                int column = 3 * j + l;
                B[i][j][3 * k + l] = dRdF[9 * column + row];
              }

        // four pointers to a 3-vector
        double * minv[4] = { &MInverse[el][0], &MInverse[el][4], &MInverse[el][8], &MInverse[el][12] }; // the four rows of MInverse (last column ignored)

        double dRdx[108]; // derivative of the element rotation matrix with respect to the positions of the tet vertices; column-major
        for(int k=0; k<4; k++)
          for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
            {
              double temp[3];
              MATRIX_VECTOR_MULTIPLY3X3(B[i][j], minv[k], temp);
              int row = 3 * i;
              int column = 3 * k + j;
              VECTOR_SET3(&dRdx[9 * column + row], temp);
            }

        // add contribution of dRdx to KElement

        // term 1: \hat{dR/dxl} K (R^T x - m)

        // compute K (R^T x - m)
        double tempVec[12]; // R^T x - m
        for(int vtx=0; vtx<4; vtx++)
        {
          double pos[3];
          for(int i=0; i<3; i++)
            pos[i] = P[4 * i + vtx];
          MATRIX_VECTOR_MULTIPLY3X3T(R, pos, &tempVec[3*vtx]);
          // subtract m
          for(int i=0; i<3; i++)
            tempVec[3*vtx+i] -= undeformedPositions[3 * vtxIndex[vtx] + i];
        }
        double a[12]; // a = K * tempVec
        for (int i=0; i<12; i++)
        {
          a[i] = 0.0;
          for (int j=0; j<12; j++)
            a[i] += KElementUndeformed[el][12 * i + j] * tempVec[j];
        }

        // add [\hat{dR/dxl} K R^T x]_l, l=1 to 12
        for(int column=0; column<12; column++)
        {
          double b[12]; // b = \hat{dR/dxl} * a
          for(int j=0; j<4; j++)
          {
            MATRIX_VECTOR_MULTIPLY3X3(&dRdx[9 * column], &a[3*j], &b[3*j]);
          }
          // write b into KElement (add b to i-th column)
          for(int row=0; row<12; row++)
            KElement[12 * row + column] += b[row]; // KElement is row-major
        }

        // term 2: (R K \hat{dRdxl}^T)x

        // re-write positions into a
        for(int vtx=0; vtx<4; vtx++)
        {
          for(int i=0; i<3; i++)
            a[3 * vtx + i] = P[4 * i + vtx];
        }

        // compute [\hat{dRdxl}^T x)]_l, l=1 to 12
        for(int column=0; column<12; column++)
        {
          double b[12]; // b = \hat{dRdxl}^T * a
          for(int j=0; j<4; j++)
          {
            MATRIX_VECTOR_MULTIPLY3X3T(&dRdx[9 * column], &a[3*j], &b[3*j]);
          }

          // add RK * b to column of KElement
          int rowStart = 0;
          for (int row=0; row<12; row++)
          {
            double contrib = 0.0;
            for (int j=0; j<12; j++)
              contrib += RK[rowStart + j] * b[j];
            KElement[rowStart + column] += contrib;
            rowStart += 12;
          }
        }
      }
    }
    else
    {
      // no warp
      memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);
      // f = K u
      double fElement[12];
      for(int i=0; i<12; i++)
      {
        fElement[i] = 0;
        for(int j=0; j<4; j++)
        {
          fElement[i] += 
            KElement[12 * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
            KElement[12 * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
            KElement[12 * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
        }
      }

      // add fElement into the global f
      if (f != NULL)
      {
        for(int j=0; j<4; j++)
        {
          f[3 * vtxIndex[j] + 0] += fElement[3 * j + 0];
          f[3 * vtxIndex[j] + 1] += fElement[3 * j + 1];
          f[3 * vtxIndex[j] + 2] += fElement[3 * j + 2];
        }
      }
    }

    if (stiffnessMatrix != NULL)
    {
      int * rowIndex = rowIndices[el];
      int * columnIndex = columnIndices[el];

      // add KElement to the global stiffness matrix
      for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
          for(int k=0; k<3; k++)
            for(int l=0; l<3; l++)
              stiffnessMatrix->AddEntry(3 * rowIndex[i] + k, 3 * columnIndex[4 * i + j] + l, KElement[12 * (3 * i + k) + 3 * j + l]);
    }
  }
}

void CorotationalLinearFEM::ClearRowColumnIndices()
{
  for (int el=0; el < tetMesh->getNumElements(); el++)
  {
    free(rowIndices[el]);
    free(columnIndices[el]);
  }

  free(rowIndices);
  free(columnIndices);
}

void CorotationalLinearFEM::BuildRowColumnIndices(SparseMatrix * sparseMatrix)
{
  int numElements = tetMesh->getNumElements();

  rowIndices = (int**) malloc (sizeof(int*) * numElements);
  columnIndices = (int**) malloc (sizeof(int*) * numElements);

  for (int el=0; el < numElements; el++)
  {
    // the 4 rows corresponding to the 4 vertices
    rowIndices[el] = (int*) malloc (sizeof(int) * 4);
    for(int i=0; i<4; i++)
      rowIndices[el][i] = tetMesh->getVertexIndex(el, i);

    // the 4 columns corresponding to all 4 vertices, in row of each vertex
    columnIndices[el] = (int*) malloc (sizeof(int) * 16);
    // find index of vertex j in row of vertex i, and cache it
    for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        columnIndices[el][4 * i + j] = sparseMatrix->GetInverseIndex(3*rowIndices[el][i], 3*rowIndices[el][j]) / 3;
  }
}

// inverse of a 3x3 matrix
// row-major format
void CorotationalLinearFEM::inverse3x3(double * A, double * AInv)
{
  // converted to C from Mathematica output   
  AInv[0] = -A[5] * A[7] + A[4] * A[8]; 
  AInv[1] = A[2] * A[7] - A[1] * A[8]; 
  AInv[2] = -A[2] * A[4] + A[1] * A[5];
  AInv[3] = A[5] * A[6] - A[3] * A[8]; 
  AInv[4] = -A[2] * A[6] + A[0] * A[8]; 
  AInv[5] = A[2] * A[3] - A[0] * A[5];
  AInv[6] = -A[4] * A[6] + A[3] * A[7]; 
  AInv[7] = A[1] * A[6] - A[0] * A[7];
  AInv[8] = -A[1] * A[3] + A[0] * A[4];

  double invDet = 1.0 / (A[0] * AInv[0] + A[1] * AInv[3] + A[2] * AInv[6]);

  for(int i=0; i<9; i++)
    AInv[i] *= invDet;
}

// inverse of a 4x4 matrix
// row-major format
void CorotationalLinearFEM::inverse4x4(double * A, double * AInv)
{
  // converted to C from Mathematica output   
  AInv[0] = -A[11] * A[14] * A[5] + A[10] * A[15] * A[5] + A[11] * A[13] * A[6] - A[10] * A[13] * A[7] - A[15] * A[6] * A[9] + A[14] * A[7] * A[9];
  AInv[1] = A[1] * A[11] * A[14] - A[1] * A[10] * A[15] - A[11] * A[13] * A[2] + A[10] * A[13] * A[3] + A[15] * A[2] * A[9] - A[14] * A[3] * A[9];
  AInv[2] = -A[15] * A[2] * A[5] + A[14] * A[3] * A[5] + A[1] * A[15] * A[6] - A[13] * A[3] * A[6] - A[1] * A[14] * A[7] + A[13] * A[2] * A[7];
  AInv[3] = A[11] * A[2] * A[5] - A[10] * A[3] * A[5] - A[1] * A[11] * A[6] + A[1] * A[10] * A[7] + A[3] * A[6] * A[9] - A[2] * A[7] * A[9];
  AInv[4] = A[11] * A[14] * A[4] - A[10] * A[15] * A[4] - A[11] * A[12] * A[6] + A[10] * A[12] * A[7] + A[15] * A[6] * A[8] - A[14] * A[7] * A[8];
  AInv[5] = -A[0] * A[11] * A[14] + A[0] * A[10] * A[15] + A[11] * A[12] * A[2] - A[10] * A[12] * A[3] - A[15] * A[2] * A[8] + A[14] * A[3] * A[8];
  AInv[6] = A[15] * A[2] * A[4] - A[14] * A[3] * A[4] - A[0] * A[15] * A[6] + A[12] * A[3] * A[6] + A[0] * A[14] * A[7] - A[12] * A[2] * A[7];
  AInv[7] = -A[11] * A[2] * A[4] + A[10] * A[3] * A[4] + A[0] * A[11] * A[6] - A[0] * A[10] * A[7] - A[3] * A[6] * A[8] + A[2] * A[7] * A[8];
  AInv[8] = -A[11] * A[13] * A[4] + A[11] * A[12] * A[5] - A[15] * A[5] * A[8] + A[13] * A[7] * A[8] + A[15] * A[4] * A[9] - A[12] * A[7] * A[9];
  AInv[9] = -A[1] * A[11] * A[12] + A[0] * A[11] * A[13] + A[1] * A[15] * A[8] - A[13] * A[3] * A[8] - A[0] * A[15] * A[9] + A[12] * A[3] * A[9];
  AInv[10] = -A[1] * A[15] * A[4] + A[13] * A[3] * A[4] + A[0] * A[15] * A[5] - A[12] * A[3] * A[5] + A[1] * A[12] * A[7] - A[0] * A[13] * A[7];
  AInv[11] = A[1] * A[11] * A[4] - A[0] * A[11] * A[5] + A[3] * A[5] * A[8] - A[1] * A[7] * A[8] - A[3] * A[4] * A[9] + A[0] * A[7] * A[9]; 
  AInv[12] = A[10] * A[13] * A[4] - A[10] * A[12] * A[5] + A[14] * A[5] * A[8] - A[13] * A[6] * A[8] - A[14] * A[4] * A[9] + A[12] * A[6] * A[9];
  AInv[13] = A[1] * A[10] * A[12] - A[0] * A[10] * A[13] - A[1] * A[14] * A[8] + A[13] * A[2] * A[8] + A[0] * A[14] * A[9] - A[12] * A[2] * A[9]; 
  AInv[14] = A[1] * A[14] * A[4] - A[13] * A[2] * A[4] - A[0] * A[14] * A[5] + A[12] * A[2] * A[5] - A[1] * A[12] * A[6] + A[0] * A[13] * A[6];
  AInv[15] = -A[1] * A[10] * A[4] + A[0] * A[10] * A[5] - A[2] * A[5] * A[8] + A[1] * A[6] * A[8] + A[2] * A[4] *A[9] - A[0] * A[6] * A[9];

  double invDet = 1.0 / (A[0] * AInv[0] + A[1] * AInv[4] + A[2] * AInv[8] + A[3] * AInv[12]);

  for(int i=0; i<16; i++)
    AInv[i] *= invDet;
}

