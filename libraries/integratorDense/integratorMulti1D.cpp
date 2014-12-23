/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "integrator" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC     *
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

#include <math.h>
#include <memory.h>
#include "matrix.h"
#include "matrixLAPACK.h"
#include "integratorMulti1D.h"
#include "float.h"

IntegratorMulti1D::IntegratorMulti1D(int r, double timestep, double * massMatrix_, double * tangentStiffnessMatrix_, double dampingMassCoef, double dampingStiffnessCoef): IntegratorBaseDense(r, timestep, dampingMassCoef, dampingStiffnessCoef)
{
  z = (double *) calloc (r, sizeof(double));
  zvel = (double *) calloc (r, sizeof(double));
  zaccel = (double *) calloc (r, sizeof(double));
  
  double * lambda = (double*) malloc (sizeof(double) * r);
  Q = (double *) malloc (sizeof(double) * r * r);

  // must do deep copy since the matrix may be modified inside the SymmetricMatrixGeneralEigenDecomposition
  double * tempTangentStiffnessMatrix = (double *) malloc (sizeof(double) * r * r);
  double * tempMassMatrix = (double *) malloc (sizeof(double) * r * r);
  memcpy(tempTangentStiffnessMatrix, tangentStiffnessMatrix_, sizeof(double) * r * r); 
  memcpy(tempMassMatrix, massMatrix_, sizeof(double) * r * r);

  SymmetricMatrixGeneralEigenDecomposition(r, tempTangentStiffnessMatrix, tempMassMatrix, Q, lambda);  // rotation is stored in Q after this function call

  // set the tangent stiffness matrix (tangent stiffness matrix is a diagonal matrix; therefore only the diagonal elements are stored.)
  tangentStiffnessMatrix = (double *) malloc (sizeof(double) * r);
  for(int dim=0; dim<r; dim++)
    tangentStiffnessMatrix[dim] = lambda[dim];
  free(lambda);

  tangentStiffnessMatrixOriginal = (double *) malloc (sizeof(double) * r);
  memcpy(tangentStiffnessMatrixOriginal, tangentStiffnessMatrix, sizeof(double) * r);

  // mass matrix is not needed because Q^T M Q = I
  massMatrix = NULL;
  
  dampingMatrix = (double *) malloc (sizeof(double) * r);

  // compute Q^T M
  QTM = NULL;
  ComputeQTM(massMatrix_);
 
  // testing: Q^T * M * Q = I 
  /*
  printf("testing: Q^T * M * Q = I\n");
  Matrix<double> QTMM(r, r, QTM, false, false);
  Matrix<double> QM(r, r, Q, false, false);  
  double * ITemp = (double *) malloc (sizeof(double) * r * r);
  Matrix<double> ITempM(r, r, ITemp, false, false);
  ITempM = QTMM * QM;
  printf("\nQ^T * M * Q = \n");
  for(int row=0; row<r; row++)
  {
    for(int col=0; col<r; col++)
    {
      double temp = ITemp[ELT(r, row, col)];
      if (temp < 1e-8)
        printf("0, ");
      else
        printf("%G, ", temp);
    }
    printf("\n");
  }
  free(ITemp);
  */

  rotatedForces = (double *) malloc (sizeof(double) * r);

  free(tempTangentStiffnessMatrix);
  free(tempMassMatrix);
}

IntegratorMulti1D::IntegratorMulti1D(int r, double timestep, double * massMatrix_, double * tangentStiffnessMatrixDiagonalElements, double * modeRotationMatrix, double dampingMassCoef, double dampingStiffnessCoef)
: IntegratorBaseDense(r, timestep, dampingMassCoef, dampingStiffnessCoef)
{
  z = (double *) calloc (r, sizeof(double));
  zvel = (double *) calloc (r, sizeof(double));
  zaccel = (double *) calloc (r, sizeof(double));

  Q = (double *) malloc (sizeof(double) * r * r);
  memcpy(Q, modeRotationMatrix, sizeof(double) * r * r);
  
  tangentStiffnessMatrix = (double *) malloc (sizeof(double) * r);
  memcpy(tangentStiffnessMatrix, tangentStiffnessMatrixDiagonalElements, sizeof(double) * r);

  tangentStiffnessMatrixOriginal = (double *) malloc (sizeof(double) * r);
  memcpy(tangentStiffnessMatrixOriginal, tangentStiffnessMatrix, sizeof(double) * r);

  massMatrix = NULL;
  dampingMatrix = (double *) malloc (sizeof(double) * r);

  // compute (Q^T)M
  QTM = NULL;
  ComputeQTM(massMatrix_);

  rotatedForces = (double *) malloc (sizeof(double) * r);
}

IntegratorMulti1D::~IntegratorMulti1D()
{
  free(z);
  free(zvel);
  free(zaccel);
  free(tangentStiffnessMatrixOriginal);
  free(Q);
  free(QTM);
  free(rotatedForces);
}

void IntegratorMulti1D::ResetToRest()
{
  IntegratorBase::ResetToRest();
  memset(z, 0, sizeof(double) * r);
  memset(zvel, 0, sizeof(double) * r);
  memset(zaccel, 0, sizeof(double) * r);
}

void IntegratorMulti1D::RotateVector(double * Q, double * source, double * dest, bool useRotationTranspose)
{  
  Matrix<double> sourceM(r, 1, source, false, false);
  double * result = (double *) malloc (sizeof(double) * r);
  Matrix<double> resultM(r, 1, result, false, false);
  Matrix<double> QM(r, r, Q, false, false);
  
  if (useRotationTranspose)
    resultM = QM.MultiplyT(sourceM);
  else
    resultM = QM * sourceM;

  memcpy(dest, result, sizeof(double) * r);
  free(result);
}

int IntegratorMulti1D::SetState(double * q_, double * qvel_)
{
  memcpy(q, q_, sizeof(double) * r);
  RotateVector(QTM, q, z);

  if (qvel_ == NULL)
  {
    memset(qvel, 0, sizeof(double) * r);
    memset(zvel, 0, sizeof(double) * r);
  }
  else
  {
    memcpy(qvel, qvel_, sizeof(double) * r);
    RotateVector(QTM, qvel, zvel);
  }

  // zaccel + (dampingMassCoef + dampingStiffnessCoef * tangentStiffnessMatrixDiagonalElement) * zvel + tangentStiffnessMatrixDiagonalElement * z = Q^T * f_ext
  // rotate the external forces: rotatedForces = Q^T * externalForces
  bool useRotationTranspose = true;
    RotateVector(Q, externalForces, rotatedForces, useRotationTranspose);

  for(int dim=0; dim<r; dim++)
  {
    double zvelCoef = dampingMassCoef + dampingStiffnessCoef * tangentStiffnessMatrix[dim];
    zaccel[dim] = rotatedForces[dim] - tangentStiffnessMatrix[dim] * z[dim] - zvelCoef * zvel[dim];
  }
  return 0;
}

void IntegratorMulti1D::ConstrainToSphere(double R2)
{
  double norm2 = 0;
  for(int i=0; i<r; i++)
    norm2 += q[i] * q[i];

  if (norm2 > R2)
  {
    double beta = sqrt(R2 / norm2);
    for(int dim=0; dim<r; dim++)
    {
      q[dim] *= beta;
      q_1[dim] = q[dim];
      z[dim] *= beta;

      qvel[dim] = 0;
      qvel_1[dim] = 0;
      zvel[dim] = 0;

      qaccel[dim] = 0;
      qaccel_1[dim] = 0;
      zaccel[dim] = 0;
    }
  }
}

void IntegratorMulti1D::SetqState(const double * q_, const double * qvel_, const double * qaccel_)
{
  // Since q = Qz and (Q^T)MQ = I, z = (Q^T)Mq

  if (q != q_)
    memcpy(q, q_, sizeof(double) * r);
  RotateVector(QTM, q, z);
  
  if (qvel_ != NULL)
  {
    if (qvel != qvel_)
      memcpy(qvel, qvel_, sizeof(double) * r);
    RotateVector(QTM, qvel, zvel);
  }

  if (qaccel_ != NULL)
  {
    if (qaccel != qaccel_)
      memcpy(qaccel, qaccel_, sizeof(double) * r); 
    RotateVector(QTM, qaccel, zaccel);
  }
}

void IntegratorMulti1D::GetqState(double * q_, double * qvel_, double * qaccel_)
{
  // q = Qz
  RotateVector(Q, z, q);
  memcpy(q_, q, sizeof(double) * r);
 
  if (qvel_ != NULL)
  {
    RotateVector(Q, zvel, qvel);
    memcpy(qvel_, qvel, sizeof(double) * r);
  }
  if (qaccel_ != NULL) 
  {
    RotateVector(Q, zaccel, qaccel);
    memcpy(qaccel_, qaccel, sizeof(double) * r);
  }
}

void IntegratorMulti1D::SetQ(int index, double qIndex)
{
  q[index] = qIndex;
  RotateVector(QTM, q, z);
}

double IntegratorMulti1D::GetQ(int index)
{
  RotateVector(Q, z, q);
  return q[index];
}

double * IntegratorMulti1D::Getq()
{
  RotateVector(Q, z, q);
  return q;
}

double * IntegratorMulti1D::Getqvel()
{
  RotateVector(Q, zvel, qvel);
  return qvel;
}

double * IntegratorMulti1D::Getqaccel()
{
  RotateVector(Q, zaccel, qaccel);
  return qaccel;
}

void IntegratorMulti1D::ComputeQTM(double * massMatrix_)
{
  free(QTM);
  QTM = (double *) malloc (sizeof(double) * r * r);
  Matrix<double> QM(r, r, Q, false, false);
  Matrix<double> QTMM(r, r, QTM, false, false);
  Matrix<double> massM(r, r, massMatrix_, false, false);
  QTMM = QM.MultiplyT(massM);
}

