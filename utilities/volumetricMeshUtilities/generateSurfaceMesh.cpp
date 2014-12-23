/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "generateSurfaceMesh" utility , Copyright (C) 2007 CMU, 2009 MIT,     *
 *                                               2014 USC                *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this utility in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "getopts.h"
#include "volumetricMeshLoader.h"
#include "generateSurfaceMesh.h"

/*
  Creates the surface mesh of the given volumetric mesh.
  Note: interior volumetric mesh vertices are kept in the surface mesh (as isolated vertices).
        So, the vertex set of the volumetric mesh is identical to the surface mesh vertex set,
        with the same order.
*/

int main( int argc, char** argv )
{

  if ( argc < 3 )
  {
    std::cout << "Usage: " << argv[0] << " [volumetric mesh file] [output obj file] [-t]" << std::endl;
    std::cout << "Generates the surface mesh of the given volumetric mesh." << std::endl;
    std::cout << "-t: (this option only applies with cubic meshes) Triangulates each quad face. Default: faces are left as quads." << std::endl;
    return 1;
  }

  char * meshFile = argv[1];
  char * outputFile = argv[2];

  bool outputTriangleMesh = false;

  opt_t opttable[] =
  {
    { (char*)"t", OPTBOOL, &outputTriangleMesh },
    { NULL, 0, NULL }
  };

  argv += 1;
  argc -= 1;
  getopts(argc,argv,opttable);

  if (outputTriangleMesh)
    printf("Triangle mesh will be triangulated (if needed).\n");

  VolumetricMesh * mesh = VolumetricMeshLoader::load(meshFile);

  GenerateSurfaceMesh generateSurfaceMesh;
  ObjMesh * objMesh = generateSurfaceMesh.ComputeMesh(mesh, outputTriangleMesh);
  objMesh->save(outputFile);

  return 0;
}

