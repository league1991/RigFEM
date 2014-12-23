/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2014 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
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

#ifndef MASSSPRINGSYSTEMFROMCUBICMESH_H_
#define MASSSPRINGSYSTEMFROMCUBICMESH_H_

#include "massSpringSystem.h"

/*
  This class can generate a mass-spring system from the given cube mesh: 
  each cube mesh vertex becomes a mass point, and each pair of vertices
  in a cube is connected by an edge.
  See also massSpringSystem.h
*/

class CubicMesh;

class MassSpringSystemFromCubicMesh
{
public:

  static int GenerateMassSpringSystem(CubicMesh * CubicMesh, MassSpringSystem ** massSpringSystem, double density, double tensileStiffness, double damping, int addGravity=0);
};

#endif

