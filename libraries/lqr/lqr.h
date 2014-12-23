/*************************************************************************
 *                                                                       *
 * "lqr" library , Copyright (C) 2009 MIT                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Jovan Popovic                                *
 * Funding: Singapore-MIT GAMBIT Game Lab                                *
 * Version 1.0                                                           *
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

#ifndef _LQR_H_
#define _LQR_H_

/*
  The "LQR" class implements a linear-quadratic regulator (LQR)
  for the following linear time-varying dynamical system:

  \dot{x} = F(t) * x + G(t) * u   (1),

  where x is the state, u is the control, and F and G are (potentially time-varying) matrices. 
  The state and the control can have arbitrary dimensions (not necessarily equal).

  The system in equation (1) is continuous: both x and u are viewed as continuous signals in time.
  In practice, the system is discretized in time, leading to a sequence of states x_0, x_1, ..., x_T,
  and control vecotrs u_0, u_1, ..., u_{T-1}. Here, T is the total number of timesteps.
  The "LQR" class computes optimal control for the **continuous** system of (1),
  under the following assumptions:
  1. The applied control is constant over each timestep.
  2. The matrices F(t) and G(t) are constant over each timestep 

  In practice, computing LQR under these assumptions requires evaluating certain definite integrals
  involving matrix exponential functions, for which no analytical formulas exist. This code uses
  Simpson integration to approximate these integrals (see the function "BuildRiccatiMatrices").

  LQR is perhaps the simplest kind of an optimal controller for a dynamical system.
  When the dynamical system is of the form in Equation (1), the LQR controller is optimal,
  in the following sense. At any timestep i, it produces control that minimizes the following energy functional:
  E = 1/2 * integral from i*h to T*h of x(tau)^T * Q(tau) * x(tau) +
      1/2 * x(T*h)^T * Q_f * x(T*h) +
      1/2 * integral from i*h to T*h of u(tau)^T * R(tau) * u(tau) ,
  where h is the timestep, and Q, Q_f and R are the position error cost, final position error cost, and control cost, respectively. 

  The system (1) is first-order. Second-order dynamical systems
  can be transcribed into the first-order format by augmenting
  the state to (x, \dot{x}) (analogous for higher-order).

  The LQR controller can also be used to control nonlinear systems,
  by building a controller around a trajectory of a nonlinear system.
  One example can be seen in [1].

  [1] Jernej Barbic, Jovan Popovic: Real-time Control of Physically Based Simulations using Gentle Forces, ACM Transactions on Graphics 27(5) (SIGGRAPH Asia 2008), Singapore, Dec 2008
  [2] STENGEL, R. F. 1994. Optimal Control and Estimation. Dover Publications, New York.
*/

#include <vector>
#include "matrix.h"

class LQR
{
public:
  
  // builds a LQR controller
  // Fs and Gs are the matrices F_i and G_i (see above), for i=0,1,...T-1.
  // cost matrices are of the form lambda * I
  LQR(int stateSize, int controlSize, double dt,
      double positionCost, double controlCost, double terminalCost,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs,
      int verbose=1);

  // same as above except the cost matrices are general (dense)
  // note: cost matrices are assumed to be time-independent (code could be modified to support time-varying cost matrices without too much difficulty)
  LQR(int stateSize, int controlSize, double dt,
      double * positionCost, double * controlCost, double * terminalCost,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs, int verbose=1);

  // load controller from file
  LQR(char * filename); 

  ~LQR();

  // save controller to file
  void Save(char * filename);

  // compute the optimal control "du", at timestep "timestepIndex", when the system is in state "dx"
  void ComputeControl(int timestepIndex, double * dx, double * du);

  // compute the objective function value, if the system is in state "dx" at timestep 0
  double ComputeInitialCost(double * dx);

  // returns the gain matrix, i.e., the matrix K that relates:
  // du = K * dx
  Matrix<double> & GetGainMatrix(int timestepIndex);

protected:
  // the size of the state space
  int stateSize;
  // the size of the control input
  int controlSize;
  // the time step
  double dt;

  // feedback matrices
  std::vector<Matrix<double> > Ks; // du = K * dx

  // store value matrix at t=0
  Matrix<double> * P0;

  // builds the controller
  void Init(int stateSize, int controlSize, double dt,
      double * positionCost, double * controlCost, double * terminalCost,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs, int verbose=1);

  void BuildRiccatiMatrices(double dt, double * F, double * G,
    Matrix<double> & Q, Matrix<double> & R,
    Matrix<double> & Phik, Matrix<double> & Gammak,
    Matrix<double> & Qk, Matrix<double> & Mk, Matrix<double> & Rk);
    
  bool Nan_check();
  bool my_isnan(double x);
};

#endif

