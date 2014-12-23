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

#ifndef _INTEGRATORBASE_H_
#define _INTEGRATORBASE_H_

/*

This code can numerically timestep a system of ODEs of the form:

M * q'' + (alpha * M + beta * K(q)) q' + R(q) = fext(t) , where

q is an unknown vector function, M is an arbitrary symmetric positive-definite 
matrix, fext(t) are arbitrary user-provided external forces, alpha and beta 
are arbitrary non-negative scalar constants (Rayleigh damping parameters), 
R(q) is an arbitrary user-provided vector function, and K(q) = d R / d q is 
its user-provided gradient. 

Such a system arises, for example, when simulating a nonlinear
deformable object using the Finite Element Method (FEM), or also using
mass-spring systems.  This code has been used for simulations of large 
deformations of 3D solid deformable objects, with dynamics (see [1]).

The code supports several numerical integrators (see derived classes).

The code can handle both large sparse systems and dense reduced systems.
For example, the large sparse version can be used to simulate a general 
deforming tetrahedral mesh. The dense version can be used for simulations 
that employ model reduction (see [1]).

The class in this file (IntegratorBase) is the abstract base class. In practice, 
you need to use one of the provided derived classes, depending on whether 
you want to simulate large sparse systems or dense reduced systems.

All these classes are generic in that you can provide your own arbitrary 
internal forces R(q) and their gradients K(q). You do so by deriving from the 
InternalForceModel class, and passing that class to the appropriate integrator 
class constructor. Several force model classes are provided, including 
the cubic polynomial reduced StVK model from [1], and a linearized version of 
that model.

For dense simulations, you need a BLAS and LAPACK library. This is necessary 
to solve the r x r linear systems inside the implicit Newmark solver.
We have successfully used the following BLAS and LAPACK libraries:
1. Windows: Intel Math Kernel Library for Windows
2. Red Hat Linux: Intel Math Kernel Library for Linux
3. Mac OS X: both BLAS and LAPACK are already included with Mac OS X 
   (the vecLib framework)
You can probably also download the generic BLAS and LAPACK implementations 
from www.netlib.org. All dense matrices are stored in column-major format.

For large sparse simulations, you need a large sparse linear system solver. 
The code supports the following solvers:
1. SPOOLES (a free solver), 
2. PARDISO (we used the commercial version that comes with Intel MKL; this 
solver is multi-threaded and can be executed across multiple cores of a CPU),
3. our own Jacobi-preconditioned Conjugate Solver (from our "sparseMatrix" library).
We were able to run sparse simulations on Windows, Linux and Mac OS X.

The code also supports static simulations, i.e., simulations where the dynamic
terms are neglected, and the system only computes the static equilibrium
under the currently applied external forces.

References:
[1] Jernej Barbic, Doug L. James: Real-Time Subspace Integration for 
St.Venant-Kirchhoff Deformable Models, ACM Transactions on Graphics 24(3) 
(SIGGRAPH 2005), p. 982-990, Los Angeles, CA, August 2005
[2] Jernej Barbic: Real-time Reduced Large-Deformation Models and Distributed 
Contact for Computer Graphics and Haptics, PhD Thesis, Carnegie Mellon University, 
August 2007

Both publications are available online at www.jernejbarbic.com .

*/

#include <stdlib.h>

// This abstract class is derived into: IntegratorBaseDense (dense systems)
// and ImplicitNewmarkSparse ((large) sparse systems).
class IntegratorBase
{
public:
  // r is the dimension of the simulation; it equals 3*n for unreduced systems, where n is the number of vertices in the simulation mesh; with reduction, r equals the size of the simulation basis 
  // the damping coefficients are tangential Rayleigh damping coefficients, see [2]
  IntegratorBase(int r, double timestep, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0);

  virtual ~IntegratorBase();

  // === set/get the state (position,velocity and acceleration) ===

  // set integrator to zero (q, qvel, and qaccel are set to zero)
  virtual void ResetToRest(); 

  // sets the position and the velocity
  // this routine will internally automatically compute proper acceleration 
  // returns 0 on success, 1 if solver fails to converge
  virtual int SetState(double * q, double * qvel=NULL) = 0;

  // sets the position, velocity, and acceleration
  // note: if you don't set all three at once, old values will persist, which may not be what you want
  virtual void SetqState(const double * q, const double * qvel=NULL, const double * qaccel=NULL);

  // copies the state into spaces provided by q,qvel,qaccel (each a vector of length r; if NULL is provided for either of q,qvel,qccel, that part of the state is not copied)
  virtual void GetqState(double * q, double * qvel=NULL, double * qaccel=NULL);

  // set/get invidivual position components:
  inline virtual void SetQ(int index, double qIndex) { q[index] = qIndex; } 
  inline virtual double GetQ(int index) { return q[index]; } 

  // obtain pointers to the internally stored position, velocity, and acceleration
  // (advanced usage)
  inline virtual double * Getq() { return q; }
  inline virtual double * Getqvel() { return qvel; }
  inline virtual double * Getqaccel() { return qaccel; }

  // == set external forces (a vector of r numbers) ===

  // external forces remain in force until explicity changed
  void SetExternalForces(double * externalForces); 
  void AddExternalForces(double * externalForces); 
  void GetExternalForces(double * externalForces);
  double * GetExternalForces() { return externalForces; };
  void SetExternalForcesToZero(); 

  // === set integration and simulation parameters ===

  inline virtual void SetTimestep(double timestep) { this->timestep = timestep; }
  inline double GetTimeStep() { return timestep; }

  // scale all internal forces and stiffness matrix entries by internalForceScalingFactor (e.g., to make the model softer or stiffer; default: 1.0)
  // note: frequency spectrum is scaled by sqrt(internalForceScalingFactor)
  virtual void SetInternalForceScalingFactor(double internalForceScalingFactor) { this->internalForceScalingFactor = internalForceScalingFactor; }

  // tangential Rayleigh damping parameters
  inline virtual void SetDampingMassCoef(double dampingMassCoef) { this->dampingMassCoef = dampingMassCoef; }
  inline virtual void SetDampingStiffnessCoef(double dampingStiffnessCoef) { this->dampingStiffnessCoef = dampingStiffnessCoef;}
  inline double GetDampingMassCoef() { return dampingMassCoef; }
  inline double GetDampingStiffnessCoef() { return dampingStiffnessCoef; }

  // === perform one timestep of the simulation ===
  // (the choice of integrator depends on the derived class)

  virtual int DoTimestep() = 0;

  // === misc ===

  inline int GetNumDOFs() { return r; }
  inline int Getr() { return r; }

  virtual double GetKineticEnergy() = 0;
  virtual double GetTotalMass() = 0; // sum of mass matrix entries

  virtual double GetForceAssemblyTime() = 0;
  virtual double GetSystemSolveTime() = 0;

  // constrain the system to ||q||^2 < R2
  // useful to prevent large values from occuring
  virtual void ConstrainToSphere(double R2);

protected:

  double * q; // current deformation amplitudes
  double * qvel; // current velocities of deformation amplitudes
  double * qaccel; // current acceleration (used inside implicit newmark integration)
  double * qresidual, * qdelta; // aux integration variables 
  double * q_1; // deformation amplitudes at previous time-step
  double * qvel_1;
  double * qaccel_1;

  double * internalForces; // current internal force amplitudes
  double * externalForces; // current external force amplitudes

  double internalForceScalingFactor; 

  double * buffer;

  // these two store the damping parameters 
  double dampingMassCoef;
  double dampingStiffnessCoef;

  int r; // number of reduced DOFs 

  double timestep; 
};

#endif

