/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * Merges several obj files into one.                                    *
 * Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                            *
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
  Merges several obj files into one.
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
using namespace std;
#include "getopts.h"
#include "matrixIO.h"
#include "objMesh.h"

vector<string> filenames;

int main( int argc, char** argv )
{
  const int numFixedArguments = 3;
  if( argc < numFixedArguments )
  {
    std::cout << "Merges several obj meshes into one. This routine uses the \"objMesh\" class to load the meshes. It can also merge materials." << std::endl;
    std::cout << "Usage: " << argv[0] << " [list of obj files] [output file] [-m]" << std::endl;
    std::cout << "  -m : also output materials" << std::endl;
    std::cout << "  -d : (if -m) also remove duplicate materials" << std::endl;
    return 1;
  }

  bool outputMaterials = false; 
  bool removeDuplicatedMaterials = false; 
  char * objListFilename = argv[1];
  char * objMeshname = argv[2];

  opt_t opttable[] =
  {
    { (char*)"m", OPTBOOL, &outputMaterials },
    { (char*)"d", OPTBOOL, &removeDuplicatedMaterials },
    { NULL, 0, NULL }
  };
  
  argv += numFixedArguments-1;
  argc -= numFixedArguments-1;
  int optup = getopts(argc,argv,opttable);
  if (optup != argc)
  {
    printf("Error parsing options. Error at option %s.\n",argv[optup]);
    return 1;
  }

  ObjMesh outputObjMesh;

  char s[4096];
  FILE * fin;
  char fileMode[96] = "r";
  OpenFile_(objListFilename, &fin, fileMode);

  while(fgets(s,4096,fin) != 0)
  {
    if(s[strlen(s)-1] == '\n')
      s[strlen(s)-1] = '\0';

    printf("%s\n",s);

    ObjMesh * objMesh = new ObjMesh(s);

    // add vertices
    int numVerticesCurrent = outputObjMesh.getNumVertices();
    for(unsigned int i=0; i<objMesh->getNumVertices(); i++)
      outputObjMesh.addVertexPosition(objMesh->getPosition(i));

    // add normals
    int numNormalsCurrent = outputObjMesh.getNumNormals();
    for(unsigned int i=0; i<objMesh->getNumNormals(); i++)
      outputObjMesh.addVertexNormal(objMesh->getNormal(i));
   
    // add texture coordinates
    int numTextureCoordinatesCurrent = outputObjMesh.getNumTextureCoordinates();
    for(unsigned int i=0; i<objMesh->getNumTextureCoordinates(); i++)
      outputObjMesh.addTextureCoordinate(objMesh->getTextureCoordinate(i));

    // add materials
    int numMaterialsCurrent = outputObjMesh.getNumMaterials();
    for(unsigned int i=0; i<objMesh->getNumMaterials(); i++)
      outputObjMesh.addMaterial(objMesh->getMaterial(i));

    for(unsigned int i=0; i<objMesh->getNumGroups(); i++)
    {
      const ObjMesh::Group * group = objMesh->getGroupHandle(i);
      outputObjMesh.addGroup(group->getName());
      unsigned int newGroupID = outputObjMesh.getNumGroups() - 1;
      ObjMesh::Group * newGroup = (ObjMesh::Group*) outputObjMesh.getGroupHandle(newGroupID);
      newGroup->setMaterialIndex(numMaterialsCurrent + group->getMaterialIndex());

      // over all faces in the group of the current obj file
      for(unsigned int j=0; j<group->getNumFaces(); j++)
      {
        const ObjMesh::Face * face = group->getFaceHandle(j);
        for(unsigned int k=0; k<face->getNumVertices(); k++)
        {
          ObjMesh::Vertex * vertex = (ObjMesh::Vertex*) face->getVertexHandle(k);
          vertex->setPositionIndex(vertex->getPositionIndex() + numVerticesCurrent);
          if (vertex->hasNormalIndex())
            vertex->setNormalIndex(vertex->getNormalIndex() + numNormalsCurrent);
          if (vertex->hasTextureCoordinateIndex())
            vertex->setTextureCoordinateIndex(vertex->getTextureCoordinateIndex() + numTextureCoordinatesCurrent);
        }
        outputObjMesh.addFaceToGroup(*face,newGroupID);
      }
    }

    delete(objMesh);
  }
  
  fclose(fin);

  if (outputMaterials)
  {
    if (removeDuplicatedMaterials)
    {
      printf("Removing duplicated materials..."); fflush(NULL);
      int numNewMaterials = outputObjMesh.removeDuplicatedMaterials();
      printf(" retained %d materials.\n", numNewMaterials);
    }
    outputObjMesh.save(objMeshname, 1);
  }
  else
  {
    outputObjMesh.save(objMeshname);
  }

  return(0);
}

