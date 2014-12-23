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

#include "tetMesh.h"
#include "massSpringSystemFromTetMesh.h"

int MassSpringSystemFromTetMesh::GenerateMassSpringSystem(TetMesh * tetMesh, MassSpringSystem ** massSpringSystem, double density, double tensileStiffness, double damping, int addGravity)
{
  int numParticles;
  double * restPositions;
  int numTets;
  int * tets;
  int numVerticesPerElement;

  tetMesh->exportMeshGeometry(&numParticles, &restPositions, &numTets, &numVerticesPerElement, &tets);

  if (numVerticesPerElement != 4)
  {
    printf("Sanity check error: mesh is not a tet mesh.\n");
    return 1;
  }

  //*massSpringSystem = new MassSpringSystem(numParticles, restPositions, numTets, tets, density, tensileStiffness, damping, addGravity);
  *massSpringSystem = new MassSpringSystem(numParticles, restPositions, TET, numTets, tets, density, tensileStiffness, damping, addGravity);

  free(tets);
  free(restPositions);

  return 0;
}

