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
  A class to timestep dense diagonal dynamics.
*/

#ifndef _IMPLICITNEWMARKDENSEMULTI1D_H_
#define _IMPLICITNEWMARKDENSEMULTI1D_H_

#include "integratorMulti1D.h"

class ImplicitNewmarkDenseMulti1D : public IntegratorMulti1D, public ImplicitNewmarkDense
{
public:
  ImplicitNewmarkDenseMulti1D(int r, double timestep, double * massMatrix, double * tangentStiffnessMatrix, double dampingMassCoef=0.0, double dampingStiffnessCoef=0.0, double NewmarkBeta=0.25, double NewmarkGamma=0.5);

  virtual ~ImplicitNewmarkDenseMulti1D();
  
  virtual void SetTimestep(double timestep);
  virtual void SetNewmarkBeta(double NewmarkBeta);
  virtual void SetNewmarkGamma(double NewmarkGamma);
  virtual void SetInternalForceScalingFactor(double internalForceScalingFactor);
  virtual void SetDampingMassCoef(double dampingMassCoef);
  virtual void SetDampingStiffnessCoef(double dampingStiffnessCoef);

  // performs one timestep of simulation (returns 0 on success, and 1 on failure)
  virtual int DoTimestep(); 

protected:
  double * coef_deltaZ;     // a r x 1 vector, for dimension i, coef_deltaZ[i] = alpha1 + diagonalizedDampingMatrix[i] * alpha4 + diagonalizedStiffness[i]
  double * coef_zvel;       // a r x 1 vector, for dimension i, coef_zvel[i] = alpha2 - diagonalizedDampingMatrix[i] * alpha5
  double * coef_zaccel;     // a r x 1 vector, for dimension i, coef_zaccel[i] = alpha3 - diagonalizedDampingMatrix[i] * alpha6
  
  void UpdateCoefs();       // this function should be always called after alphas are updated
};

#endif

