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

#ifndef _ISOTROPICHYPERELASTICFEMMT_H_
#define _ISOTROPICHYPERELASTICFEMMT_H_

#include "isotropicHyperelasticFEM.h"

/*
  This class is a multi-threaded version of the class "IsotropicHyperelasticFEM".
  It uses POSIX threads ("pthreads") as the threading API.
  Each thread assembles the internal force with respect to a subset of all the mesh elements. 
  At the end, the individual results are added into a global internal force vector.

  See also "isotropicHyperelasticFEM.h".
*/

class IsotropicHyperelasticFEMMT : public IsotropicHyperelasticFEM
{
public:
  // see "isotropicHyperelasticFEM.h" for usage
  // numThreads is the number of threads to use for the computation
  IsotropicHyperelasticFEMMT(TetMesh * tetMesh, IsotropicMaterial * isotropicMaterial, double principalStretchThreshold=-DBL_MAX, bool addGravity=false, double g=9.81, int numThreads=1);
  virtual ~IsotropicHyperelasticFEMMT();

  // Computes strain energy, internal forces, and/or tangent stiffness matrix, as requested by computationMode. It returns 0 on success, and non-zero on failure.
  // this is an advanced function; you normally do not need to use it
  virtual int GetEnergyAndForceAndTangentStiffnessMatrixHelper(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode);

  int GetStartElement(int rank);
  int GetEndElement(int rank);

protected:
  int numThreads;
  int * startElement, * endElement;
  double * energyBuffer;
  double * internalForceBuffer;
  SparseMatrix ** tangentStiffnessMatrixBuffer;

  void Initialize();
};

#endif

