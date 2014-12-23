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

#include "objMesh.h"
#include "configFile.h"
#include "massSpringSystemFromObjMesh.h"
#include "massSpringSystemFromObjMeshConfigFile.h"
#include <string.h>

char * MassSpringSystemFromObjMeshConfigFile::DuplicateString(const char * s)
{
  // strdup sometimes causes problems, so we use this
  char * p = (char*) malloc (sizeof(char) * (strlen(s) + 1));
  memcpy(p, s, sizeof(char) * (strlen(s) + 1));
  return p;
}

int MassSpringSystemFromObjMeshConfigFile::GenerateMassSpringSystem(const char * configFilename, MassSpringSystem ** massSpringSystem, MassSpringSystemObjMeshConfiguration * massSpringSystemObjConfiguration)
{
  char massSpringMeshFilename[4096];
  double surfaceDensity, tensileStiffness, shearStiffness, bendStiffness, damping;
  int addGravity;

  printf("Parsing configuration file %s...\n", configFilename);
  ConfigFile configFile;
  configFile.addOption("massSpringMeshFilename", massSpringMeshFilename);
  configFile.addOption("surfaceDensity", &surfaceDensity);
  configFile.addOption("tensileStiffness", &tensileStiffness);
  configFile.addOption("shearStiffness", &shearStiffness);
  configFile.addOption("bendStiffness", &bendStiffness);
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

  ObjMesh * quadMesh;
  try
  {
    quadMesh = new ObjMesh(massSpringMeshFilename);
  }
  catch(int eCode)
  {
    printf("Error: unable to load mesh from %s. Code: %d\n", massSpringMeshFilename, eCode);
    return 1;
  };

  MassSpringSystemFromObjMesh massSpringSystemFromObjMesh;

  int code = massSpringSystemFromObjMesh.GenerateMassSpringSystem(quadMesh, massSpringSystem, surfaceDensity, tensileStiffness, shearStiffness, bendStiffness, damping, addGravity);
  delete(quadMesh);

  if (massSpringSystemObjConfiguration != NULL)
  {
    massSpringSystemObjConfiguration->massSpringMeshFilename = DuplicateString(massSpringMeshFilename);
    massSpringSystemObjConfiguration->surfaceDensity = surfaceDensity;
    massSpringSystemObjConfiguration->tensileStiffness = tensileStiffness;
    massSpringSystemObjConfiguration->shearStiffness = shearStiffness;
    massSpringSystemObjConfiguration->bendStiffness = bendStiffness;
    massSpringSystemObjConfiguration->damping = damping;
    massSpringSystemObjConfiguration->addGravity = addGravity;
  }

  return code;
}

