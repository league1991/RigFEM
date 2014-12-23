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
#include "configFile.h"
#include "massSpringSystemFromTetMesh.h"
#include "massSpringSystemFromTetMeshConfigFile.h"

char * MassSpringSystemFromTetMeshConfigFile::DuplicateString(const char * s)
{
  // strdup sometimes didn't work well, so we used this
  char * p = (char*) malloc (sizeof(char) * (strlen(s) + 1));
  memcpy(p, s, sizeof(char) * (strlen(s) + 1));
  return p;
}

int MassSpringSystemFromTetMeshConfigFile::GenerateMassSpringSystem(const char * configFilename, MassSpringSystem ** massSpringSystem, MassSpringSystemTetMeshConfiguration * massSpringSystemTetMeshConfiguration)
{
  char tetMeshFilename[4096];
  char surfaceMeshFilename[4096];
  double density, tensileStiffness, damping;
  int addGravity;

  printf("Parsing configuration file %s...\n", configFilename);
  ConfigFile configFile;
  configFile.addOption("tetMeshFilename", tetMeshFilename);
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

  TetMesh * tetMesh;
  try
  {
    tetMesh = new TetMesh(tetMeshFilename);
  }
  catch(int eCode)
  {
    printf("Error: unable to load mesh from %s. Code: %d\n", tetMeshFilename, eCode);
    return 1;
  };
  printf("Tet mesh loaded.\n");

  MassSpringSystemFromTetMesh massSpringSystemFromTetMesh;

  int code = massSpringSystemFromTetMesh.GenerateMassSpringSystem(tetMesh, massSpringSystem, density, tensileStiffness, damping, addGravity);

  delete(tetMesh);

  if (massSpringSystemTetMeshConfiguration != NULL)
  {
    massSpringSystemTetMeshConfiguration->tetMeshFilename = DuplicateString(tetMeshFilename);
    massSpringSystemTetMeshConfiguration->surfaceMeshFilename = DuplicateString(surfaceMeshFilename);
    massSpringSystemTetMeshConfiguration->density = density;
    massSpringSystemTetMeshConfiguration->tensileStiffness = tensileStiffness;
    massSpringSystemTetMeshConfiguration->damping = damping;
    massSpringSystemTetMeshConfiguration->addGravity = addGravity;
  }

  return code;
}

