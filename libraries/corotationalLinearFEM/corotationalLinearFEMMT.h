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

#ifndef _COROTATIONALLINEARFEMMT_H_
#define _COROTATIONALLINEARFEMMT_H_

#include "corotationalLinearFEM.h"

/*
   Multi-threaded version of the CorotationalLinearFEM class. 
   It uses the POSIX threads ("pthreads").
   See also corotationalLinearFEM.h
*/

class CorotationalLinearFEMMT : public CorotationalLinearFEM
{
public:

  CorotationalLinearFEMMT(TetMesh * tetMesh, int numThreads=1);
  virtual ~CorotationalLinearFEMMT();

  virtual void ComputeForceAndStiffnessMatrix(double * vertexDisplacements, double * internalForces, SparseMatrix * stiffnessMatrix, int warp=1);

  int GetStartElement(int rank);
  int GetEndElement(int rank);

protected:
  int numThreads;
  int * startElement, * endElement;
  double * internalForceBuffer;
  SparseMatrix ** stiffnessMatrixBuffer;

  void Initialize();
  void ComputeHelper(double * u, double * uSecondary, void * target, bool addQuantity);

};

#endif

