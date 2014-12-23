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

#include "cubicMesh.h"
#include "massSpringSystemFromCubicMesh.h"

int MassSpringSystemFromCubicMesh::GenerateMassSpringSystem(CubicMesh * cubicMesh, MassSpringSystem ** massSpringSystem, double density, double tensileStiffness, double damping, int addGravity)
{
  int numParticles;
  double * restPositions;
  int numCubes;
  int * cubes;
  int numVerticesPerElement;

  cubicMesh->exportMeshGeometry(&numParticles, &restPositions, &numCubes, &numVerticesPerElement, &cubes);

  if (numVerticesPerElement != 8)
  {
    printf("Sanity check error: mesh is not a cube mesh.\n");
    return 1;
  }

  *massSpringSystem = new MassSpringSystem(numParticles, restPositions, CUBE, numCubes, cubes, density, tensileStiffness, damping, addGravity);

  free(cubes);
  free(restPositions);

  return 0;
}

