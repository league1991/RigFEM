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
#include "configFile.h"
#include "massSpringSystemFromCubicMesh.h"
#include "massSpringSystemFromCubicMeshConfigFile.h"

char * MassSpringSystemFromCubicMeshConfigFile::DuplicateString(const char * s)
{
  // strdup sometimes didn't work well, so we used this
  char * p = (char*) malloc (sizeof(char) * (strlen(s) + 1));
  memcpy(p, s, sizeof(char) * (strlen(s) + 1));
  return p;
}

int MassSpringSystemFromCubicMeshConfigFile::GenerateMassSpringSystem(const char * configFilename, MassSpringSystem ** massSpringSystem, MassSpringSystemCubicMeshConfiguration * massSpringSystemCubicMeshConfiguration)
{
  char cubicMeshFilename[4096];
  char surfaceMeshFilename[4096];
  double density, tensileStiffness, damping;
  int addGravity;

  printf("Parsing configuration file %s...\n", configFilename);
  ConfigFile configFile;
  configFile.addOption("cubicMeshFilename", cubicMeshFilename);
  configFile.addOptionOptional("surfaceMeshFilename", surfaceMeshFilename, "__none");
  configFile.addOption("density", &density);
  configFile.addOption("tensileStiffness", &tensileStiffness);
  configFile.addOption("damping", &damping);
  configFile.addOption("addGravity", &addGravity);

  if (configFile.parseOptions(configFilename) != 0)
  {
    printf("Error parsing options.\n");
    return 1;
  }

  // the config variables have now been loaded with their specified values

  // informatively print the variables (with assigned values) that were just parsed
  configFile.printOptions();

  CubicMesh * cubeMesh;
  try
  {
    cubeMesh = new CubicMesh(cubicMeshFilename);
  }
  catch(int eCode)
  {
    printf("Error: unable to load mesh from %s. Code: %d\n", cubicMeshFilename, eCode);
    return 1;
  };
  printf("sq mesh loaded.\n");

  MassSpringSystemFromCubicMesh massSpringSystemFromCubicMesh;

  int code = massSpringSystemFromCubicMesh.GenerateMassSpringSystem(cubeMesh, massSpringSystem, density, tensileStiffness, damping, addGravity);

  delete(cubeMesh);

  if (massSpringSystemCubicMeshConfiguration != NULL)
  {
    massSpringSystemCubicMeshConfiguration->cubicMeshFilename = DuplicateString(cubicMeshFilename);
    massSpringSystemCubicMeshConfiguration->surfaceMeshFilename = DuplicateString(surfaceMeshFilename);
    massSpringSystemCubicMeshConfiguration->density = density;
    massSpringSystemCubicMeshConfiguration->tensileStiffness = tensileStiffness;
    massSpringSystemCubicMeshConfiguration->damping = damping;
    massSpringSystemCubicMeshConfiguration->addGravity = addGravity;
  }

  return code;
}

