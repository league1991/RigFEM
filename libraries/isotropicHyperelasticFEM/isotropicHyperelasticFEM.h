/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "isotropic hyperelastic FEM" library , Copyright (C) 2014 USC         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin                            *
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

#ifndef _ISOTROPICHYPERELASTICFEM_H_
#define _ISOTROPICHYPERELASTICFEM_H_

#include <float.h>
#include "tetMesh.h"
#include "sparseMatrix.h"
#include "isotropicMaterial.h"

/*
  Implementation of hyperelastic isotropic nonlinear FEM elasticity, using
  linear tetrahedral elements. 

  This code can compute strain energy, internal forces and tangent stiffness matric
es (derivative of internal forces). It can also handle and restore from element inversion.

  The tangent stiffness matrix is computed using the block-diagonal dP/dF tensor, as described in [Teran 2005].

  Citations:

  IRVING G., TERAN J., FEDKIW R.: Invertible Finite
  Elements for Robust Simulation of Large Deformation. In Proc.  of the Symp. on Comp. Animation 2004 (2004), pp. 131–140.

  TERAN J., SIFAKIS E., IRVING G., FEDKIW R.: Robust
  Quasistatic Finite Elements and Flesh Simulation. In 2005
  ACM SIGGRAPH / Eurographics Symp. on Computer Animation
  (July 2005), pp. 181–190.

*/

class IsotropicHyperelasticFEM
{
public:
  // Before creating this class, you must first create the tet mesh, and create an instance of the "IsotropicMaterial" material (e.g., NeoHookeanMaterial).
  // 
  // The inversionThreshold controls when the inversion prevention mechanism activates:
  // If the principal stretches (the "lambdas") are smaller than the inversionThreshold, they will be clamped to that, which will generally prevent element inversion. For example, a typical inversionThreshold value would be 0.1. By default, inversion handling is disabled (inversionThreshold=-infinity). Values of 1.0 or higher should not be used.
  //
  // The material properties are determined as follows:
  // If the "isotropicMaterial" is of the "Homogeneous" kind, then the material properties are homogeneous and are specified by "isotropicMaterial"; material properties in the "tetMesh" are ignored (only geometry is used)
  // If the "isotropicMaterial" is not of the "Homogeneous" kind, then the material properties are such as specified in the "tetMesh" and "isotropicMaterial" class (and may be non-homogeneous)
  IsotropicHyperelasticFEM(TetMesh * tetMesh, IsotropicMaterial * isotropicMaterial, double inversionThreshold=-DBL_MAX, bool addGravity=false, double g=9.81);
  virtual ~IsotropicHyperelasticFEM();

  double ComputeEnergy(double * u); // get the nonlinear elastic strain energy

  // get the nonlinear internal forces
  // both vertex displacements "u" and internal forces refer to the vertices of the simulation mesh
  // they must be (pre-allocated) vectors of length 3 * numVertices
  // the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
  // i.e., the computed internal forces are negatives of the actual physical internal forces acting on the material
  void ComputeForces(double * u, double * internalForces);

  // allocate memory for the non-zero entries of the stiffness matrix
  void GetStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix);
  // get the nonlinear stiffness matrix given the vertex displacement vector u
  void GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix);
  // get both nonlinear internal forces and nonlinear stiffness matrix
  void GetForceAndTangentStiffnessMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix);

  // compute damping forces based on the velocity of the vertices,
  // see p6 section 6.2 of [Irving 04]
  void ComputeDampingForces(double dampingPsi, double dampingAlpha, double * u, double * uvel, double * dampingForces);

  // enables or disables the gravity (note: you can also set this in the constructor; use this routine to turn the gravity on/off during the simulation)
  void SetGravity(bool addGravity) { this->addGravity = addGravity; } // if addGravity is enabled, ComputeForces will subtract the gravity force from the internal forces (note: subtraction, not addition, is used because the internal forces are returned with the sign as described in the f_int(x) comment above)

  inline TetMesh * GetTetMesh() { return tetMesh; }

  void SetMaterial(IsotropicMaterial * isotropicMaterial_) { isotropicMaterial = isotropicMaterial_; }

  // === Advanced functions below; you normally do not need to use them: ===
  // Computes strain energy, internal forces, and/or tangent stiffness matrix, as requested by computationMode. It returns 0 on success, and non-zero on failure.
  // computationMode:
  // bit 0: compute energy
  // bit 1: compute internal force
  // bit 2: compute stiffness matrix
  typedef enum { COMPUTE_ENERGY=1, COMPUTE_INTERNALFORCES=2, COMPUTE_TANGENTSTIFFNESSMATRIX=4 } computationModeType;
  virtual int GetEnergyAndForceAndTangentStiffnessMatrixHelper(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode);
  // Initialization for "GetEnergyAndForceAndTangentStiffnessMatrixPrologue" (must always be called before calling "GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse")
  void GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode);
  // The workhorse (main computational routine); processes mesh elements startEl <= el < endEl (assembles partial strain energy, internal forces, and/or tangent stiffness matrix, as requested by computationMode. It returns 0 on success, and non-zero on failure.
  int GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(int startEl, int endEl, double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode);

protected:
  TetMesh * tetMesh; // the tet mesh
  IsotropicMaterial * isotropicMaterial; // the material 

  // acceleration indices
  int ** row_;
  int ** column_;

  double * restVerticesPosition;    // length equals to the #vertices in the mesh times 3
  double * currentVerticesPosition; // it equals restVerticesPosition + u

  // If a principal stretch (i.e., the F^hat in Irving 2004 paper) is smaller than 
  // inversionThreshold, it will be clamped to that.  This is important to ensure invertibility. 
  // For example, a typical value for invertible StVK would be 0.5.
  double inversionThreshold;

  bool addGravity;
  double g;

  // this is the b=(A1N1 + A2N2 + A3N3) in the paper,
  // see p.3 section 4
  Vec3d * areaWeightedVertexNormals;
  // this is an array of inv(Dm), Dm is the matrix
  // of the edge vectors of a tet. See p3 section 3 of [Irving 04]
  // length of this array equals to the number of tet in the mesh
  Mat3d * dmInverses;
  // an array of deformation gradient F
  Mat3d * Fs;
  // an array of F^hat (i.e., the principal stretches)
  Vec3d * Fhats;
  // an array of the V rotation matrices
  Mat3d * Vs;
  // an array of the U rotation matrices
  Mat3d * Us;

  // tet volumes; necessary to compute the elastic strain energy
  double * tetVolumes;
  void ComputeTetVolumes();

  // the area weighted vertex normals in the rest configuration
  void ComputeAreaWeightedVertexNormals();
  // compute inv(Dm)
  void PrepareDeformGrad(); // called once in the constructor

  // Compute strain energy based on the user-specified material
  virtual double ComputeEnergyFromStretches(int elementIndex, double * lambda);
  // Compute the diagonalized first Piola-Kirchhoff stress P^hat
  virtual void ComputeDiagonalPFromStretches(int elementIndex, double * lambda, double * PDiag);
  // Compute the element stiffness matrix
  virtual void ComputeTetK(int el, double K[144], int clamped);
  // Compute the derivative of the first Piola Kirchhoff stress P with respect to
  // the deformation gradient F. Since P and F both have 9 entries, dPdF has 81 entries
  virtual void Compute_dPdF(int el, double dPdF[81], int clamped);
  // Compute the derivative of the deformation gradient F with respect 
  // to the displacement vector u
  void Compute_dFdU();
  // The G is a 3x3 matrix where the columns are the nodal forces 
  // (see p3 section 4 of [Irving 04]). So dGdF is the derivative of G (i.e., nodal forces)
  // with respect to the deformation gradient F
  void Compute_dGdF(Vec3d * b0, Vec3d * b1, Vec3d * b2, double dPdF[], double dGdF[]);

  // {i,j,m} goes from 0 to 2 inclusively,
  // and {n} goes from 0 to 3 inclusively.
  // converts 3x3x3x4 tensor indices to 9x12 matrix indices
  int tensor9x12Index(int i, int j, int m, int n);

  // dFdUs is an array of dFdU (i.e., derivative of the deformation gradient with
  // respect to the displacement vector u), and dFdU is stored as a array of doubles.
  double * dFdUs; // array of length 9x12 x numElements
  // Ds is the matrix which the columns are the edge vector of a tet (see p3 
  // section 3 of [Irving 04]). dDSdU is a 9x12 matrix which stores the derivative
  // of the Ds matrix with respect to the displacement vector u. Because Ds has 9 entries
  // and the u vector (of a single tet) has 12 entries (i.e., 4 vertices * 3 dof), so I 
  // re-arrange dDSdU to be a 9x12 matrix
  double dDSdU[108]; //in 9x12 matrix format

  void dP_From_dF(Mat3d & dF, Mat3d & dP);
  //this gammaValue function is used by dP_dF
  inline double gammaValue(int i, int j, double sigma[3], double invariants[3], double gradient[3], double hessian[6]);

  // {i,j,m,n} goes from 0 to 2 inclusively
  // converts 3x3x3x3 tensor indices to 9x9 row-major matrix indices
  int tensor9x9Index(int i, int j, int m, int n);

  /*
    rowMajorMatrixToTeran is a mapping from a row-major 3x3 matrix to a 9-vector in Teran's order
    i.e., 
             | m00 m01 m02 |   | 0 1 2 |
    from M = | m10 m11 m12 | = | 3 4 5 |
             | m20 m21 m22 |   | 6 7 8 |

    to V = [m00 m11 m22 m01 m10 m02 m20 m12 m21] = [0 4 8 1 3 2 6 5 7]
  */
  int rowMajorMatrixToTeran[9];

  /*
    "teranToRowMajorMatrix" is a mapping from a 9-vector (in Teran's order) to row-major 3x3 matrix

    from V = [m00 m11 m22 m01 m10 m02 m20 m12 m21]

           | m00 m01 m02 |   | 0 1 2 |
    to M = | m10 m11 m12 | = | 3 4 5 |
           | m20 m21 m22 |   | 6 7 8 |
  */
  int teranToRowMajorMatrix[9];
};

#endif

