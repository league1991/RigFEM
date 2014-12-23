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

/*
  A base abstract class to timestep decoupled dynamics
  For each dimension, the decoupled equation is as follows,
  
  z_accel + d * z_vel + k * z = f_external

  (i.e., usually linear model reduction dynamics).

  It is derived into classes that support decoupled implicit Newmark,
  and decoupled central differences.

*/

#ifndef _INTEGRATORMULTI1D_H_
#define _INTEGRATORMULTI1D_H_

#include "integratorBaseDense.h"
#include "reducedForceModel.h"

class IntegratorMulti1D : public virtual IntegratorBaseDense
{
public:
  // r is the dimension of the simulation basis
  // the damping coefficients are tangential Rayleigh damping coefficients, see [2]

  // Constructor 1: Mass matrix and tangent stiffness matrix are two general r x r matrices. The constructor will diagonalize the system making it decoupled.
  IntegratorMulti1D(int r, double timestep, double * massMatrix, double * tangentStiffnessMatrix, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0);

  // Constructor 2:
  // Besides r and timestep, you need to pass:
  // massMatrix (general mass matrix, not diagonalized, r x r) which is used to compute Q^T * M;
  // tangentStiffnessMatrixDiagonalElements (diagonalized, r x 1);
  // modeRotationMatrix (denoted as Q) is a r x r mode rotation matrix, it is usually obtained by mass-PCA, that is, to solve
  // K * ksi_i = lambda_i * M * ksi_i; (i = 1, 2, ... r), ksi_i is the i-th eigen vector whose dimension is r.
  // Q = [ksi_1, ksi_2, ..., ksi_r];
  // Also in this case make sure Q^T * M * Q = I
  IntegratorMulti1D(int r, double timestep, double * massMatrix, double * tangentStiffnessMatrixDiagonalElements, double * modeRotationMatrix, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0); 

  virtual ~IntegratorMulti1D();
 
  // sets the position and the velocity
  // this routine will internally automatically compute proper acceleration 
  // returns 0 on success, 1 if solver fails to converge
  virtual int SetState(double * q, double * qvel=NULL);
 
  // sets the position, velocity, and acceleration
  // note: if you don't set all three at once, old values will persist, which may not be what you want
  virtual void SetqState(const double * q, const double * qvel=NULL, const double * qaccel=NULL);

  // copies the state into spaces provided by q,qvel,qaccel (each a vector of length r; if NULL is provided for either of q,qvel,qccel, that part of the state is not copied)
  virtual void GetqState(double * q, double * qvel=NULL, double * qaccel=NULL);

  // set/get invidivual position components:
  virtual void SetQ(int index, double qIndex);
  virtual double GetQ(int index);

  // obtain pointers to the internally stored position, velocity, and acceleration
  // (advanced usage)
  virtual double * Getq();
  virtual double * Getqvel();
  virtual double * Getqaccel();
   
  virtual void ResetToRest();

  // projects the current state so that ||q||^2 <= R2
  virtual void ConstrainToSphere(double R2);

protected: 
  double * tangentStiffnessMatrixOriginal;
  // source and dest are r x 1 vectors;
  // Q is a r x r rotation matrix;
  // dest = Rotation * source if useRotationTranspose = false
  // dest = Rotation^T * source otherwise
  void RotateVector(double * Q, double * source, double * dest, bool useRotationTranspose = false);
  double * Q;
  double * z;
  double * zvel;
  double * zaccel;

  void ComputeQTM(double * massMatrix);
  double * QTM;             // a r x r matrix, QTM = (Q^T) * massMatrix (not in pratical use)
  double * rotatedForces;   // a r x 1 vector, rotatedForces = Q^T * externalForces  
};

#endif

