/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "clothBW" library , Copyright (C) 2014 USC                            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Andy Pierce, Yu Yu Xu, Jernej Barbic                     *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
using namespace std;
#include "clothBW.h"
#include "macros.h"
#include "matrixMultiplyMacros.h"
#include "minivector.h"
//#include "performanceCounter.h"

// constructor without triangleUVs
ClothBW::ClothBW(int numParticles_, double * masses_, double * restPositions_, int numTriangles_, int * triangles_, int * triangleGroups_, int numMaterialGroups_, double * groupTensileStiffness_, double * groupShearStiffness_, double * groupBendStiffnessU_, double * groupBendStiffnessV_, double * groupDamping_, int addGravity_) : bu(1.0), bv(1.0), addGravity(addGravity_), g(9.81)
{
  double * triangleUVs_ = (double *) malloc( sizeof(double) * 3 * 2 * numTriangles_ );
  
  //default to computing everything
  _computationConditions[0] = true; // stretch&shear force
  _computationConditions[1] = true; // bend force
  _computationConditions[2] = true; // stretch&shear bend stiffness
  _computationConditions[3] = true; // bend stiffness matrix
  
  // default to taking rest angles into account
  useRestAnglesForBendingForces = 1;
  
  for (int i = 0 ; i < numTriangles_; i++)
  {	
    int particleA = triangles_[3*i+0];
    int particleB = triangles_[3*i+1];
    int particleC = triangles_[3*i+2];
    
    triangleUVs_[6*i+0] = 0.0;		// for the vertex A u
    triangleUVs_[6*i+1] = 0.0;		// for the vertex A v
    
    Vec3d x0; // vector from A to B (previously B to A)
    x0[0] = restPositions_[3*particleA+0] - restPositions_[3*particleB+0];
    x0[1] = restPositions_[3*particleA+1] - restPositions_[3*particleB+1]; 
    x0[2] = restPositions_[3*particleA+2] - restPositions_[3*particleB+2];  
    
    double lengthAB = len(x0);
    
    triangleUVs_[6*i+2] = lengthAB;	// for the vertex B u
    triangleUVs_[6*i+3] = 0.0;		// for the vertex B v
    
    Vec3d xn0 = norm(x0); //vector from A to B normalized
    
    Vec3d x1; // vector from A to C (previously C to A)
    x1[0] = restPositions_[3*particleA+0] - restPositions_[3*particleC+0];
    x1[1] = restPositions_[3*particleA+1] - restPositions_[3*particleC+1]; 
    x1[2] = restPositions_[3*particleA+2] - restPositions_[3*particleC+2];    
    
    double lengthAC = len(x1);
    
    triangleUVs_[6*i+4] = dot(xn0, x1); // for the vertex C u
    triangleUVs_[6*i+5] = sqrt(lengthAC * lengthAC - triangleUVs_[6*i+4] * triangleUVs_[6*i+4]);	// for the vertex C v
    
  }
  GenerateBW(numParticles_, masses_, restPositions_, numTriangles_, triangles_, triangleUVs_, triangleGroups_, numMaterialGroups_, groupTensileStiffness_, groupShearStiffness_, groupBendStiffnessU_, groupBendStiffnessV_, groupDamping_, addGravity_);
}

// constructor with trianglesUVs
ClothBW::ClothBW(int numParticles_, double * masses_, double * restPositions_, int numTriangles_, int * triangles_, double * triangleUVs_, int * triangleGroups_, int numMaterialGroups_, double * groupTensileStiffness_, double * groupShearStiffness_, double * groupBendStiffnessU_, double * groupBendStiffnessV_, double * groupDamping_, int addGravity_) : bu(1.0), bv(1.0), addGravity(addGravity_), g(9.81)
{
  //default to computing everything
  _computationConditions[0] = true; // stretch&shear force
  _computationConditions[1] = true; // bend force
  _computationConditions[2] = true; // stretch&shear bend stiffness
  _computationConditions[3] = true; // bend stiffness matrix
  
  // default to taking rest angles into account
  useRestAnglesForBendingForces = 1;
  
  GenerateBW(numParticles_, masses_, restPositions_, numTriangles_, triangles_, triangleUVs_, triangleGroups_, numMaterialGroups_, groupTensileStiffness_, groupShearStiffness_, groupBendStiffnessU_, groupBendStiffnessV_, groupDamping_, addGravity_);
}

// copy constructor
ClothBW::ClothBW(ClothBW & ClothBW)
{
  _computationConditions[0] = ClothBW._computationConditions[0]; // stretch&shear force
  _computationConditions[1] = ClothBW._computationConditions[1]; // bend force
  _computationConditions[2] = ClothBW._computationConditions[2]; // stretch&shear bend stiffness
  _computationConditions[3] = ClothBW._computationConditions[3]; // bend stiffness matrix
  
  // default to taking rest angles into account
  useRestAnglesForBendingForces = ClothBW.useRestAnglesForBendingForces;
  
  numParticles = ClothBW.numParticles;
  
  masses = (double*) malloc (sizeof(double) * numParticles);
  memcpy(masses, ClothBW.masses, sizeof(double) * numParticles);
  restPositions = (double*) malloc (sizeof(double) * 3 * numParticles);
  memcpy(restPositions, ClothBW.restPositions, sizeof(double) * 3 * numParticles);
  
  numQuads = ClothBW.numQuads;
  restAngles = (double*) malloc (sizeof(double) * numQuads);
  memcpy(restAngles, ClothBW.restAngles, sizeof(double) * numQuads);
  
  numTriangles = ClothBW.numTriangles;
  triangles = (int*) malloc (sizeof(int) * 3 * numTriangles);
  memcpy(triangles, ClothBW.triangles, sizeof(int) * 3 * numTriangles);
  triangleGroups = (int*) malloc (sizeof(int) * numTriangles);
  memcpy(triangleGroups, ClothBW.triangleGroups, sizeof(int) * numTriangles);
  inverseIndicesStretchAndShear = (int*) malloc (sizeof(int) * 9 * numTriangles);
  memcpy(inverseIndicesStretchAndShear, ClothBW.inverseIndicesStretchAndShear, sizeof(int) * 9 * numTriangles);
  triangleUVs = (double *) malloc(sizeof(double) * 3 * 2 * numTriangles);
  memcpy(triangleUVs, ClothBW.triangleUVs, sizeof(double) * 3 * 2 * numTriangles);
  
  
  inverseIndicesQuad = (int*) malloc (sizeof(int) * 16 * numQuads);
  memcpy(inverseIndicesQuad, ClothBW.inverseIndicesQuad, sizeof(int) * 16 * numQuads);
  quadComponentIndices = (int*) malloc (sizeof(int) * 6 * numQuads);
  memcpy(quadComponentIndices, ClothBW.quadComponentIndices, sizeof(int) * 6 * numQuads);
  
  numMaterialGroups = ClothBW.numMaterialGroups;
  groupTensileStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupTensileStiffness, ClothBW.groupTensileStiffness, sizeof(double) * numMaterialGroups);
  groupShearStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupShearStiffness, ClothBW.groupShearStiffness, sizeof(double) * numMaterialGroups);
  groupBendStiffnessU = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupBendStiffnessU, ClothBW.groupBendStiffnessU, sizeof(double) * numMaterialGroups);
  groupBendStiffnessV = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupBendStiffnessV, ClothBW.groupBendStiffnessV, sizeof(double) * numMaterialGroups);
  groupDamping = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupDamping, ClothBW.groupDamping, sizeof(double) * numMaterialGroups);
  
  bu = ClothBW.bu;
  bv = ClothBW.bv;
  addGravity = ClothBW.addGravity;
  g = ClothBW.g;	
}

void ClothBW::GenerateBW(int numParticles_, double * masses_, double * restPositions_, int numTriangles_, int * triangles_, double * triangleUVs_, int * triangleGroups_, int numMaterialGroups_, double * groupTensileStiffness_, double * groupShearStiffness_, double * groupBendStiffnessU_, double * groupBendStiffnessV_, double * groupDamping_, int addGravity_) 
{
  numParticles = numParticles_;	
  masses = (double*) malloc (sizeof(double) * numParticles);
  memcpy(masses, masses_, sizeof(double) * numParticles);
  restPositions = (double*) malloc (sizeof(double) * 3 * numParticles);
  memcpy(restPositions, restPositions_, sizeof(double) * 3 * numParticles);
  
  numTriangles = numTriangles_;
  triangles = (int*) malloc (sizeof(int) * 3 * numTriangles);
  memcpy(triangles, triangles_, sizeof(int) * 3 * numTriangles);
  triangleUVs = (double*) malloc(sizeof(double) * 3 * 2 * numTriangles);
  memcpy(triangleUVs, triangleUVs_, sizeof(double) * 3 * 2 * numTriangles);
  triangleGroups = (int*) malloc (sizeof(int) * numTriangles);
  memcpy(triangleGroups, triangleGroups_, sizeof(int) * numTriangles);
  
  numMaterialGroups = numMaterialGroups_;
  groupTensileStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupTensileStiffness, groupTensileStiffness_, sizeof(double) * numMaterialGroups);
  groupShearStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupShearStiffness, groupShearStiffness_, sizeof(double) * numMaterialGroups);
  groupBendStiffnessU = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupBendStiffnessU, groupBendStiffnessU_, sizeof(double) * numMaterialGroups);
  groupBendStiffnessV = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupBendStiffnessV, groupBendStiffnessV_, sizeof(double) * numMaterialGroups);
  groupDamping = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupDamping, groupDamping_, sizeof(double) * numMaterialGroups);
  
  // computing bending information
  
  // with these maps we parse all edges and find which ones 'belong' to
  // two separate triangles like the one pointed out below:
  /* 
             *----------*
            / \        /
           /-> \      /
          /  -> \    /
         /     ->\  /
        *----------*
   */
  // these edges are important because they can bend
  
  // SORTED ensures that edges are paired in ascending order
  #define SORTED(i,j) ( (i) <= (j) ? make_pair((i),(j)) : make_pair((j),(i)) )
  
  // edgeInfo data structure: maps an edge (represented by a pair of vertex indices)
  // to its corresponding triangle index as well as the index of the non-edge vertex.
  // it is a vector because an edge may belong to either 1 or 2 triangles:
  // if it belongs to two triangles, the order of the index vector is as follows:
  // [index of triangle 1] [index of non-edge vertex 1] [index of triangle 2] [index of non-edge vertex 2]
  // if it belongs to only 1 triangle, then the vector contains only the first 2 entries shown above
  std::map<std::pair<int,int>, std::vector<int> > edgeInfo;
  std::map<std::pair<int,int>, std::vector<int> >::iterator iter;
  for( int i = 0 ; i < numTriangles; i++ )
  {
    int particleA = triangles_[i*3+0];
    int particleB = triangles_[i*3+1];
    int particleC = triangles_[i*3+2];
    for( int j = 0 ; j < 3; j++)
    {
      std::pair<int,int> edge;
      int nonEdgeIndex;
      if( j == 0 )	{ edge = (SORTED(particleA,particleB)); nonEdgeIndex = particleC; }
      if( j == 1 )	{ edge = (SORTED(particleB,particleC)); nonEdgeIndex = particleA; }
      if( j == 2 )	{ edge = (SORTED(particleA,particleC)); nonEdgeIndex = particleB; }
      
      // if edge is not already in the list, we insert it to end of list
      // and map it to the appropriate triangle index and non-edge vertex index
      if( (iter = edgeInfo.find( edge )) == edgeInfo.end() )
      {
        std::vector<int> indices;
        indices.push_back(i);
        indices.push_back(nonEdgeIndex);
        edgeInfo.insert( make_pair( edge, indices ));
      }
      
      // if we find it in the list (which means that edge has already been mapped
      // to a triangle index and a non-edge index), then we insert the new triangle
      // index and non-edge index to the back of the index vector
      else 
      {
        iter->second.push_back(i);
        iter->second.push_back(nonEdgeIndex);
      }
    }
  }
  
  // delete edges with only one triangles
  map<std::pair<int,int>, std::vector<int> >::iterator iter_backup;
  iter = edgeInfo.begin();
  iter_backup = iter;
  int edgeError = 0;
  while(iter != edgeInfo.end())
  {
    iter_backup++;
    
    // if edge is only mapped to one triangle and one non-edge vertex index (1+1 = 2)
    // then we delete it
    if( iter->second.size() == 2 ) 
    { 
      // edge on the border (only 1 neighboring triangle)
      edgeInfo.erase(iter); 
      iter = iter_backup; 
    }
    
    // else if edge is mapped to two triangles and two non-edge vertex indices (2+2 = 4)
    // then we advance past it
    else if( iter->second.size() == 4 )	
      iter++;
    
    // size must be either 2 or 4 because every edge should be mapped to either
    // one triangle index and one non-edge index
    // or two triangles indices and two non-edge indices (if it is an edge shared by 2 triangles)
    else	
    {
      edgeError = 1;
      printf( "Edge Error!\n" );
      break;
    }
  }
  if (edgeError != 0)
    throw 1;
  
  numQuads = edgeInfo.size();

  restAngles = (double*) malloc (sizeof(double) * numQuads);

  quadComponentIndices = (int *)malloc(sizeof(int) * 6 * numQuads);
  iter = edgeInfo.begin();
  
  // below we save the information about edges belonging to more than one
  // triangle. we store the info in the data structure quadComponentIndices, organized
  // with 6 entries for each of these edges as pictured below:
  
    /* Old ordering
            0---------5
           / \       /
         /    \  4  /
        /  2   \   /
       /        \ /
      3----------1 
   */
  
    /* New ordering
            1---------3
           / \       /
         /    \  5  /
        /  4   \   /
       /        \ /
      0----------2 
   */
  
  for(int i = 0 ; i < numQuads; i++)
  {
    // Old order below
//    quadComponentIndices[ i * 6 + 0 ] = iter->first.first;		//  index of edge vertex 1
//    quadComponentIndices[ i * 6 + 1 ] = iter->first.second;		//	index of edge vertex 2
//    quadComponentIndices[ i * 6 + 2 ] = (iter->second)[0];		//	index of triangle 1
//    quadComponentIndices[ i * 6 + 3 ] = (iter->second)[1];		//	index of triangle 1 vertex which is not on the edge
//    quadComponentIndices[ i * 6 + 4 ] = (iter->second)[2];		//	index of triangle 2
//    quadComponentIndices[ i * 6 + 5 ] = (iter->second)[3];		//	index of triangle 2 vertex which is not on the edge	
    
    // New order below
    quadComponentIndices[ i * 6 + 0 ] = (iter->second)[1];		//	index of triangle 1 vertex which is not on the edge
    quadComponentIndices[ i * 6 + 1 ] = iter->first.first;		//  index of edge vertex 1
    quadComponentIndices[ i * 6 + 2 ] = iter->first.second;		//	index of edge vertex 2
    quadComponentIndices[ i * 6 + 3 ] = (iter->second)[3];		//	index of triangle 2 vertex which is not on the edge
    quadComponentIndices[ i * 6 + 4 ] = (iter->second)[0];		//	index of triangle 1
    quadComponentIndices[ i * 6 + 5 ] = (iter->second)[2];		//	index of triangle 2
    
    iter++;
  }
  
  // build the locations of non-zero entries in the stiffness matrix
  SparseMatrixOutline skeletonOutline(numParticles);

  for (int i=0; i<numParticles; i++)
  {
    // protecting against isolated vertices
    skeletonOutline.AddEntry(i, i);
  }
  
  for(int i=0; i<numTriangles; i++)
  {
    int particleA = triangles[3*i+0];
    int particleB = triangles[3*i+1];
    int particleC = triangles[3*i+2];
    
    skeletonOutline.AddEntry(particleA, particleA);
    skeletonOutline.AddEntry(particleA, particleB);
    skeletonOutline.AddEntry(particleA, particleC);
    skeletonOutline.AddEntry(particleB, particleA);
    skeletonOutline.AddEntry(particleB, particleB);
    skeletonOutline.AddEntry(particleB, particleC);
    skeletonOutline.AddEntry(particleC, particleA);
    skeletonOutline.AddEntry(particleC, particleB);
    skeletonOutline.AddEntry(particleC, particleC);
  }

  /* Old ordering below
            A         D 
            0---------5
           / \       /
         /    \  4  /
        /  2   \   /
       /        \ /
      3----------1
      C          B  
   */
  
  /* New ordering below
            B         D 
            1---------3
           / \       /
         /    \  5  /
        /  4   \   /
       /        \ /
      0----------2
      A          C  
   */
  
  for(int i = 0 ; i < numQuads; i++)
  {

    // particles for new ordering below
    int particleA = quadComponentIndices[6*i+0];	// index of non-edge vertex for triangle 1
    int particleB = quadComponentIndices[6*i+1];	// index of vertex 1 on the edge
    int particleC = quadComponentIndices[6*i+2];	// index of vertex 2 on the edge
    int particleD = quadComponentIndices[6*i+3];	// index of non-edge vertex for triangle 2
    
    skeletonOutline.AddEntry(particleA, particleB);
    skeletonOutline.AddEntry(particleA, particleC); 
    skeletonOutline.AddEntry(particleA, particleD); 
    skeletonOutline.AddEntry(particleB, particleA);
    skeletonOutline.AddEntry(particleB, particleC); 
    skeletonOutline.AddEntry(particleB, particleD); 
    skeletonOutline.AddEntry(particleC, particleA);
    skeletonOutline.AddEntry(particleC, particleB);
    skeletonOutline.AddEntry(particleC, particleD); 
    skeletonOutline.AddEntry(particleD, particleA);
    skeletonOutline.AddEntry(particleD, particleB);
    skeletonOutline.AddEntry(particleD, particleC); 
    
    //--- computing Resting Angles
    Vec3d restPositionA(&restPositions[3*particleA]);
    Vec3d restPositionB(&restPositions[3*particleB]);
    Vec3d restPositionC(&restPositions[3*particleC]);
    Vec3d restPositionD(&restPositions[3*particleD]);
    
    Vec3d vecBA = restPositionB - restPositionA;
    Vec3d vecCA = restPositionC - restPositionA;
    Vec3d vecBD = restPositionB - restPositionD;
    Vec3d vecCD = restPositionC - restPositionD;
    Vec3d vecBC = restPositionB - restPositionC;
    
    // new math ordering below
    Vec3d NA = cross(vecCA, vecBA); // normal for the first triangle
    Vec3d NB = cross(vecBD, vecCD);	// normal for the second triangle
    Vec3d E = vecBC; // edge vector
    
    // normalized normals and edge	
    Vec3d NAn = norm(NA);
    Vec3d NBn = norm(NB);
    Vec3d En = norm(E);
    
    double cosTheta = dot(NAn, NBn);
    double sinTheta;
    
    Vec3d NANB = cross(NAn, NBn); // storing the cross product of NAn and NBn
    
    sinTheta = dot(NANB, En);
    
    restAngles[i] = atan2(sinTheta, cosTheta);
  }

  // build inverse indices for stiffness matrix access
  SparseMatrix skeleton(&skeletonOutline);
  
  inverseIndicesStretchAndShear = (int*) malloc (sizeof(int) * 9 * numTriangles);
  for(int i=0; i<numTriangles; i++)
  {
    int particleA = triangles[3*i+0];
    int particleB = triangles[3*i+1];
    int particleC = triangles[3*i+2];
    inverseIndicesStretchAndShear[9*i+0] = skeleton.GetInverseIndex(particleA, particleA);
    inverseIndicesStretchAndShear[9*i+1] = skeleton.GetInverseIndex(particleA, particleB);
    inverseIndicesStretchAndShear[9*i+2] = skeleton.GetInverseIndex(particleA, particleC);
    inverseIndicesStretchAndShear[9*i+3] = skeleton.GetInverseIndex(particleB, particleA);
    inverseIndicesStretchAndShear[9*i+4] = skeleton.GetInverseIndex(particleB, particleB);
    inverseIndicesStretchAndShear[9*i+5] = skeleton.GetInverseIndex(particleB, particleC);
    inverseIndicesStretchAndShear[9*i+6] = skeleton.GetInverseIndex(particleC, particleA);
    inverseIndicesStretchAndShear[9*i+7] = skeleton.GetInverseIndex(particleC, particleB);
    inverseIndicesStretchAndShear[9*i+8] = skeleton.GetInverseIndex(particleC, particleC);
  }		
  
  inverseIndicesQuad = (int *) malloc (sizeof(int) * 16 * numQuads );
  for(int i = 0 ; i < numQuads; i++)
  {
    
    // particles for new ordering below
    int particleA = quadComponentIndices[6*i+0];	// index of non-edge vertex for triangle 1
    int particleB = quadComponentIndices[6*i+1];	// index of vertex 1 on the edge
    int particleC = quadComponentIndices[6*i+2];	// index of vertex 2 on the edge
    int particleD = quadComponentIndices[6*i+3];	// index of non-edge vertex for triangle 2
    
    
    // old order: CABD
    inverseIndicesQuad[16*i+0] = skeleton.GetInverseIndex(particleA, particleA);
    inverseIndicesQuad[16*i+1] = skeleton.GetInverseIndex(particleA, particleB);
    inverseIndicesQuad[16*i+2] = skeleton.GetInverseIndex(particleA, particleC);
    inverseIndicesQuad[16*i+3] = skeleton.GetInverseIndex(particleA, particleD);
    
    inverseIndicesQuad[16*i+4] = skeleton.GetInverseIndex(particleB, particleA);
    inverseIndicesQuad[16*i+5] = skeleton.GetInverseIndex(particleB, particleB);
    inverseIndicesQuad[16*i+6] = skeleton.GetInverseIndex(particleB, particleC);
    inverseIndicesQuad[16*i+7] = skeleton.GetInverseIndex(particleB, particleD);
    
    inverseIndicesQuad[16*i+8] = skeleton.GetInverseIndex(particleC, particleA);
    inverseIndicesQuad[16*i+9] = skeleton.GetInverseIndex(particleC, particleB);
    inverseIndicesQuad[16*i+10] = skeleton.GetInverseIndex(particleC, particleC);
    inverseIndicesQuad[16*i+11] = skeleton.GetInverseIndex(particleC, particleD);
    
    inverseIndicesQuad[16*i+12] = skeleton.GetInverseIndex(particleD, particleA);
    inverseIndicesQuad[16*i+13] = skeleton.GetInverseIndex(particleD, particleB);
    inverseIndicesQuad[16*i+14] = skeleton.GetInverseIndex(particleD, particleC);
    inverseIndicesQuad[16*i+15] = skeleton.GetInverseIndex(particleD, particleD);
  }
}

ClothBW::~ClothBW()
{
  free(masses);
  free(restPositions);
  free(triangles);
  free(triangleGroups);
  free(groupTensileStiffness);
  free(groupShearStiffness);
  free(groupBendStiffnessU);
  free(groupBendStiffnessV);
  free(inverseIndicesStretchAndShear);
  free(inverseIndicesQuad);
  free(quadComponentIndices);
  free(groupDamping);	
}

void ClothBW::GenerateMassMatrix(SparseMatrix **M, int expanded)
{
  SparseMatrixOutline outline(expanded * numParticles);
  for(int i=0; i<numParticles; i++)
    for(int j=0; j<expanded; j++)
      outline.AddEntry(expanded*i+j, expanded*i+j, masses[i]); 
  *M = new SparseMatrix(&outline);	
}

double ClothBW::GetTriangleSurfaceArea(double *p0, double *p1, double *p2)
{
  Vec3d s0(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
  Vec3d s1(p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]);
  return 0.5 * len(cross(s0, s1));
}

void ClothBW::SetComputationMode(bool conditions[4])
{
  _computationConditions[0] = conditions[0];
  _computationConditions[1] = conditions[1];
  _computationConditions[2] = conditions[2];
  _computationConditions[3] = conditions[3];
}

void ClothBW::ComputeGravity(double * f, bool addGravity)
{
  if (!addGravity)
    memset(f, 0, sizeof(double) * 3 * numParticles);
  
  for(int i=0; i<numParticles; i++)
    f[3*i+1] += -g * masses[i]; // *** MODIFICATION : changed gravity from
                                //                    acting in -z direction to -y
                                //                    i.e. from [3*i+2] -> [3*i+1]
}

void ClothBW::ComputeForce(double *u, double *f, bool addForce)
{
  if (!addForce)	
    memset(f, 0, sizeof(double) * 3 * numParticles);
  
  AddForce(u, f, 0, numTriangles, 0, numQuads);
  
  if (addGravity)
    ComputeGravity(f, true);	
}

void ClothBW::AddForce(const double * u, double * f, int startTriangle, int endTriangle, int startQuad, int endQuad)
{
  
  if (_computationConditions[0])
    AddStretchAndShearForce(u, f, startTriangle, endTriangle);
  
  if (_computationConditions[1])
    AddBendForce(u, f, startQuad, endQuad);
  
}

void ClothBW::AddStretchAndShearForce(const double * u, double * f, int startTriangle, int endTriangle)
{
  // basic variables for all the three forces
  double Cu, Cv;					// for stretch C(x) matrix	1*1
  double Cs;						// for shear C(x) matrix	1*1
  double phiCstr[3][2][3];		// for stretch dC(x)/dxi, 3 row, 2 column, 3 numbers
  double phiCsh[3][3];			// for shear dC(x)/dxi, 3 row, 1 column, 3 numbers
  
  for(int i = startTriangle; i < endTriangle; i++)
  {
    int group = triangleGroups[i];			// group index
    
    int particleA = triangles[3*i+0];		// triangle indices, clockwise as A, B, C
    int particleB = triangles[3*i+1];
    int particleC = triangles[3*i+2];
    
    //--- compute C(x) for Stretch and Shear
    
    Vec3d positionA(&restPositions[3*particleA]);
    positionA += Vec3d(&u[3*particleA]);
    Vec3d positionB(&restPositions[3*particleB]);
    positionB += Vec3d(&u[3*particleB]);
    Vec3d positionC(&restPositions[3*particleC]);
    positionC += Vec3d(&u[3*particleC]);
    
    Vec3d vecBA = positionB - positionA;
    Vec3d vecCA = positionC - positionA;  
    
    double delta_u1, delta_u2, delta_v1, delta_v2;	// distance of neighboring particles in planar coordinates. (delta_u1, delta_v1): planar vector from A to B; (delta_u2, delta_v2): planar vector from B to A.
    delta_u1 = triangleUVs[6*i+2] - triangleUVs[6*i+0];
    delta_v1 = triangleUVs[6*i+3] - triangleUVs[6*i+1]; 
    delta_u2 = triangleUVs[6*i+4] - triangleUVs[6*i+0];
    delta_v2 = triangleUVs[6*i+5] - triangleUVs[6*i+1]; 
    
    double scale_w = 1.0/(delta_u1*delta_v2-delta_u2*delta_v1);	 // this is a intermediate scale variable while computing wu and wv
    
    Vec3d wu = (delta_v2 * vecBA - delta_v1 * vecCA) * scale_w;
    Vec3d wv = ( -delta_u2 * vecBA + delta_u1 * vecCA) * scale_w;
    
    Vec3d wun = norm(wu); //wu normalized
    Vec3d wvn = norm(wv); // wv normalized
    double length_wu = len(wu);
    double length_wv = len(wv);
    double area;		// area for specific triangle
    
    area = 0.5 * len(cross(vecBA, vecCA));
    
    Cu = area * (length_wu - bu);
    Cv = area * (length_wv - bv);
    Cs = area * dot(wu, wv);
    
    //--- compute dC(x)/dx for Stretch and Shear
    
    Vec3d dwudx;
    Vec3d dwvdx;
    dwudx[0] = (delta_v1 - delta_v2) * scale_w;
    dwudx[1] = delta_v2 * scale_w;
    dwudx[2] = -delta_v1 * scale_w;
    dwvdx[0] = (delta_u2 - delta_u1) * scale_w;
    dwvdx[1] = -delta_u2 * scale_w;
    dwvdx[2] = delta_u1 * scale_w;
    
    for(int j = 0 ; j < 3; j++)		// per number for dC(x)/dxi vector (Stretch and Shear)
      for(int m = 0 ; m < 3; m++)	// per particle in a triangle
      {
        phiCstr[m][0][j] = wun[j] * area * dwudx[m];
        phiCstr[m][1][j] = wvn[j] * area * dwvdx[m];
        phiCsh[m][j] = area * ( dwudx[m] * wv[j] + dwvdx[m] * wu[j] );
      }
    
    //--- f = - k * C(x) * dC(x)/dx
    
    double force[3];					// storing the force per particle
    for(int j = 0 ; j < 3; j++)			// loop of 3 particles inside a triangle
      for(int m = 0 ; m < 3; m++)		// loop of 3 numbers inside a 3*1 vector
      {
        force[m] = - groupTensileStiffness[group] * ( phiCstr[j][0][m] * Cu + phiCstr[j][1][m] * Cv );
        force[m] += (- groupShearStiffness[group]) * phiCsh[j][m] * Cs;
        
        // must use "-=" because the internal force is on the left-hand side of equation
        if(j == 0)	
          f[3*particleA+m] -= force[m];
        if(j == 1)	
          f[3*particleB+m] -= force[m];
        if(j == 2)	
          f[3*particleC+m] -= force[m];
      }
  }
}

void ClothBW::AddBendForce(const double * u, double * f, int startQuad, int endQuad)
{
  for(int i = startQuad ; i < endQuad; i++)
  {
    // new ordering below
    int particleA = quadComponentIndices[6*i+0];	// index of non-edge vertex for triangle 1
    int particleB = quadComponentIndices[6*i+1];	// index of vertex 1 on the edge
    int particleC = quadComponentIndices[6*i+2];	// index of vertex 2 on the edge
    int particleD = quadComponentIndices[6*i+3];	// index of non-edge vertex for triangle 2
    
    int group1 = triangleGroups[quadComponentIndices[i*6+4]];			// group for first triangle
    int group2 = triangleGroups[quadComponentIndices[i*6+5]];			// group for second triangle
    
    double kBendU = (groupBendStiffnessU[group1] + groupBendStiffnessU[group2]) * 0.5;
    double kBendV = (groupBendStiffnessV[group1] + groupBendStiffnessV[group2]) * 0.5;
    double kBend = ( kBendU + kBendV ) * 0.5;
    
    // compute C(x)
    
    // first define some positions and accompanying vectors
    // triangle pictured below
    
    /* 
            1---------3
           / \       /
         /    \  5  /
        /  4   \   /
       /        \ /
      0----------2 
   */
    
    
    Vec3d positionA(&restPositions[3*particleA]);
    positionA += Vec3d(&u[3*particleA]);
    Vec3d positionB(&restPositions[3*particleB]);
    positionB += Vec3d(&u[3*particleB]);
    Vec3d positionC(&restPositions[3*particleC]);
    positionC += Vec3d(&u[3*particleC]);
    Vec3d positionD(&restPositions[3*particleD]);
    positionD += Vec3d(&u[3*particleD]);
    
    Vec3d vecBA = positionB - positionA;
    Vec3d vecCA = positionC - positionA;
    Vec3d vecBD = positionB - positionD;
    Vec3d vecCD = positionC - positionD;
    Vec3d vecBC = positionB - positionC;
    
    Vec3d vecCB = positionC - positionB;
    Vec3d vecAC = positionA - positionC;
    Vec3d vecDB = positionD - positionB;
    
    // now find normals for triangles 1 and 2
    Vec3d n1 = cross(vecCA, vecBA);
    Vec3d n2 = cross(vecBD, vecCD);
    Vec3d edge = vecBC;
    
    // now normalize the... normals
    Vec3d n1N = norm(n1);
    Vec3d n2N = norm(n2);
    Vec3d edgeNorm = norm(edge);
    
    // calculate sin theta and cos theta
    // sin theta = (norm1 x norm2) * edgeVector BC
    // cos theta = norm1 * norm2
    double sinTheta = dot(cross(n1N, n2N), edgeNorm);
    double cosTheta = dot(n1N, n2N);
    
    // calculate C(bend)
    double cbend = atan2(sinTheta, cosTheta);
    
    // account for resting angles?
    if (useRestAnglesForBendingForces)
      cbend -= restAngles[i];
      
//    // test this
//    while (cbend < -PI)
//      cbend += PI;
//    while (cbend > PI)
//      cbend -= PI;
//    // end test this
    
    double lengthN1 = len(n1);
    double lengthN2 = len(n2);
    double lengthE = len(edge);
    
    // calculate dn^a/dx, dn^b/dx, and dn^e/dx
    // dna/dx = Ss(q^a) from Pritchard's formula
    // dnb/dx = Ss(q^b)
    // dne/dx = I_s(q^e)
    
    // defined in Pritchard's paper
    // q^a = [BC][CA][BA][0]
    // q^b = [0] [CD][BD][BC]
    // q^e = [0] [1] [-1] [0]
    
    Vec3d qa[4];
    qa[0] = vecCB;
    qa[1] = vecAC;
    qa[2] = vecBA;
    qa[3] = 0.0;
    
    Vec3d qb[4];
    qb[0] = 0.0;
    qb[1] = vecCD;
    qb[2] = vecDB;
    qb[3] = vecBC;
    
    double dn1[4][3][3];
    double dn2[4][3][3];
    double dne[4][3][3];
    
    // each derivative conceptually structured as array of 3x3 matrices: 
    // [3x3 matrix] [3x3 matrix] [3x3 matrix] [3x3 matrix]
    for (int entry = 0; entry < 4; entry++)
    {
      //constructSkewSymmetric(qa[entry], dn1[entry], 1.0/lengthN1);
      //constructSkewSymmetric(qb[entry], dn2[entry], 1.0/lengthN2);
      
      double * dn1p = &(dn1[entry][0][0]);
      SKEW_MATRIX(qa[entry], dn1p);
      double invLengthN1 = 1.0 / lengthN1;
      for(int ii=0; ii<3; ii++)
        for(int jj=0; jj<3; jj++)
          dn1[entry][ii][jj] *= invLengthN1;
      
      double * dn2p = &(dn2[entry][0][0]);
      SKEW_MATRIX(qb[entry], dn2p);
      double invLengthN2 = 1.0 / lengthN2;
      for(int ii=0; ii<3; ii++)
        for(int jj=0; jj<3; jj++)
          dn2[entry][ii][jj] *= invLengthN2;
    }
    for (int outer=0; outer<3; outer++)
      for (int inner=0; inner<3; inner++)
      {
        dne[0][outer][inner] = 0.0;
        dne[3][outer][inner] = 0.0;
        
        if (outer == inner)
        {
          dne[1][outer][inner] = 1.0/lengthE;
          dne[2][outer][inner] = -1.0/lengthE;
        }
        else
        {
          dne[1][outer][inner] = 0.0;
          dne[2][outer][inner] = 0.0;
        }
      }
    
    // calculate dSin/dxm_s, dCos/dxm_s, and f = -k * C(x) * dC(x)/dx
    Vec3d dSin[4], dCos[4];
    for (int particle=0; particle<4; particle++)
      for (int xyz=0; xyz<3; xyz++)
      {
        dCos[particle][xyz] = dot(dn1[particle][xyz], n2N) + dot(n1N, dn2[particle][xyz]);
        
        dSin[particle][xyz] = dot ( cross(dn1[particle][xyz], n2N) + cross(n1N, dn2[particle][xyz]), edgeNorm );
        dSin[particle][xyz] += dot ( cross(n1N, n2N), dne[particle][xyz] );
        
        double force = -kBend * cbend * (cosTheta * dSin[particle][xyz] - sinTheta * dCos[particle][xyz]);
        
        // must use "-=" because the internal force is on the left-hand side of equation
        if( particle == 0 )         
          f[3*particleA+xyz] -= force;
        else if( particle == 1 )	
          f[3*particleB+xyz] -= force;
        else if( particle == 2 )	
          f[3*particleC+xyz] -= force;
        else                        
          f[3*particleD+xyz] -= force;
      }
  }
}

void ClothBW::ComputeDampingForce(double *u, double *uvel, double *f, bool addForce)
{
  // unimplemented
}

void ClothBW::AddDampingForce(double *uvel, double *f, int startTriangle, int endTriangle)
{
  // unimplemented
}

void ClothBS_Sfunction( const Vec3d & in, Vec3d & out, int row, double scale )
{
  if( row == 0 )
  {
    out[0] = 0.0;
    out[1] = -in[2];
    out[2] = in[1];
  }
  if( row == 1 )
  {
    out[0] = in[2];
    out[1] = 0.0;
    out[2] = -in[0];
  }
  if( row == 2 )
  {
    out[0] = -in[1];
    out[1] = in[0];
    out[2] = 0.0;
  }
  out[0] *= scale;
  out[1] *= scale;
  out[2] *= scale;
}

void ClothBW::GenerateStiffnessMatrixTopology(SparseMatrix **K)
{
  SparseMatrixOutline KOutline(3*numParticles);
  
  for(int vtx=0; vtx < numParticles; vtx++)
  {
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        KOutline.AddEntry(3 * vtx + j, 3 * vtx + k);
  }
  
  // for per triangle
  for(int i=0; i<numTriangles; i++)
  {
    int particleA = triangles[3*i+0];
    int particleB = triangles[3*i+1];
    int particleC = triangles[3*i+2];
    
    if ((particleA < 0) || (particleA >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
    if ((particleB < 0) || (particleB >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
	if ((particleC < 0) || (particleC >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        KOutline.AddEntry(3 * particleA + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleA + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleA + j, 3 * particleC + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleC + k);
        KOutline.AddEntry(3 * particleC + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleC + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleC + j, 3 * particleC + k);
      }
  }
  
  // for per quad, in addition to per triangle topology
  for( int i = 0 ; i < numQuads; i++ )
  {
    // new ordering below
    int particleA = quadComponentIndices[6*i+0];	// index of non-edge vertex for triangle 1
    int particleB = quadComponentIndices[6*i+1];	// index of vertex 1 on the edge
    int particleC = quadComponentIndices[6*i+2];	// index of vertex 2 on the edge
    int particleD = quadComponentIndices[6*i+3];	// index of non-edge vertex for triangle 2
    
    if ((particleA < 0) || (particleA >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
    if ((particleB < 0) || (particleB >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
	if ((particleC < 0) || (particleC >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
	if ((particleD < 0) || (particleD >= numParticles))
    {
      printf("Particle error.\n");
      exit(1);
    }
    
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        // new order below
        KOutline.AddEntry(3 * particleA + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleA + j, 3 * particleC + k);
        KOutline.AddEntry(3 * particleA + j, 3 * particleD + k);
        KOutline.AddEntry(3 * particleD + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleD + j, 3 * particleC + k);
        KOutline.AddEntry(3 * particleD + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleC + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleD + k);
        KOutline.AddEntry(3 * particleC + j, 3 * particleD + k);
        
        // old order below
//        KOutline.AddEntry(3 * particleC + j, 3 * particleA + k);
//        KOutline.AddEntry(3 * particleC + j, 3 * particleB + k);
//        KOutline.AddEntry(3 * particleC + j, 3 * particleD + k);
//        KOutline.AddEntry(3 * particleD + j, 3 * particleA + k);
//        KOutline.AddEntry(3 * particleD + j, 3 * particleB + k);
//        KOutline.AddEntry(3 * particleD + j, 3 * particleC + k);
//        KOutline.AddEntry(3 * particleA + j, 3 * particleC + k);
//        KOutline.AddEntry(3 * particleB + j, 3 * particleC + k);
//        KOutline.AddEntry(3 * particleA + j, 3 * particleD + k);
//        KOutline.AddEntry(3 * particleB + j, 3 * particleD + k);
      }
  }
  
  *K = new SparseMatrix(&KOutline);
}

void ClothBW::ComputeStiffnessMatrix(double *u, SparseMatrix *K, bool addMatrix)
{
  if (!addMatrix)	
    K->ResetToZero();
  
  // in this single-threaded version, compute all triangles and all quads in 1 run
  AddStiffnessMatrix(u, K, 0, numTriangles, 0, numQuads);
}

void ClothBW::AddStiffnessMatrix(double *u, SparseMatrix *K, int startTriangle, int endTriangle,
                                 int startQuad, int endQuad)
{
  if (_computationConditions[2])
    AddStretchAndShearStiffnessMatrix(u, K, startTriangle, endTriangle);
  
  if (_computationConditions[3])
    AddBendStiffnessMatrix(u, K, startQuad, endQuad);
  
}

void ClothBW::AddStretchAndShearStiffnessMatrix(double *u, SparseMatrix *K, int startTriangle, int endTriangle)
{
  // basic variables for all the three forces
  double Cu, Cv;				// for stretch C(x) matrix
  double Cs;					// for shear C(x) matrix
  double phiCstr[3][2][3];		// for stretch dC(x)/dxi, 3 row, 2 column, 3 numbers
  double phiCsh[3][3];			// for shear dC(x)/dxi, 3 row, 1 column, 3 numbers
  
  for(int i = startTriangle; i < endTriangle; i++)
  {
    int group = triangleGroups[i];			// group index
    int particleA = triangles[3*i+0];		// triangle indices, clockwise as A, B, C
    int particleB = triangles[3*i+1];
    int particleC = triangles[3*i+2];
	
    //--- compute C(x) for Stretch and Shear
    
    Vec3d positionA(&restPositions[3*particleA]);
    positionA += Vec3d(&u[3*particleA]);
    Vec3d positionB(&restPositions[3*particleB]);
    positionB += Vec3d(&u[3*particleB]);
    Vec3d positionC(&restPositions[3*particleC]);
    positionC += Vec3d(&u[3*particleC]);
    
    Vec3d vecBA = positionB - positionA;
    Vec3d vecCA = positionC - positionA;  
    
    double delta_u1, delta_u2, delta_v1, delta_v2;	// distance of neighboring particles in planar coordinates. (delta_u1, delta_v1): planar vector from A to B; (delta_u2, delta_v2): planar vector from B to A.
    delta_u1 = triangleUVs[6*i+2] - triangleUVs[6*i+0];
    delta_v1 = triangleUVs[6*i+3] - triangleUVs[6*i+1]; 
    delta_u2 = triangleUVs[6*i+4] - triangleUVs[6*i+0];
    delta_v2 = triangleUVs[6*i+5] - triangleUVs[6*i+1]; 
    
    double scale_w = 1.0/(delta_u1*delta_v2-delta_u2*delta_v1);	 // this is a intermediate scale variable while computing wu and wv
    
    Vec3d wu = (delta_v2 * vecBA - delta_v1 * vecCA) * scale_w; 
    Vec3d wv = ( -delta_u2 * vecBA + delta_u1 * vecCA) * scale_w; 
    
    Vec3d wun = norm(wu); //wu normalized
    Vec3d wvn = norm(wv); // wv normalized
    double length_wu = len(wu);
    double length_wv = len(wv);
    double area;		// area for specific triangle
    
    area = 0.5 * len(cross(vecBA, vecCA));
    
    // preventing unstable situation
    //if(area > 0.5)	
      //printf("Warning: large area\n" );
    
    Cu = area * (length_wu - bu);
    Cv = area * (length_wv - bv);
    Cs = area * dot(wu, wv);
    
    //--- compute dC(x)/dx for Stretch and Shear
    
    Vec3d dwudx;	// dwu/dx0_x, dwu/dx0_y, dwu/dx0_z
    Vec3d dwvdx;	// dwv/dx0_x, dwv/dx0_y, dwv/dx0_z
    dwudx[0] = (delta_v1 - delta_v2) * scale_w;
    dwudx[1] = delta_v2 * scale_w;
    dwudx[2] = -delta_v1 * scale_w;
    dwvdx[0] = (delta_u2 - delta_u1) * scale_w;
    dwvdx[1] = -delta_u2 * scale_w;
    dwvdx[2] = delta_u1 * scale_w;
    
    for(int j = 0 ; j < 3; j++)		// per number for dC(x)/dxi Stretch
      for(int m = 0 ; m < 3; m++)	// per particle in a triangle
      {
        phiCstr[m][0][j] = wun[j] * area * dwudx[m]; 
        phiCstr[m][1][j] = wvn[j] * area * dwvdx[m]; 
        phiCsh[m][j] = area * ( dwudx[m] * wv[j] + dwvdx[m] * wu[j] );
      }
    
    //--- compute Stiffness Matrix for Stretch and Shear
    // Kij = df(xi) / dxj = - k * ( (dC(x) / dxi) * transpose(dC(x) / dxj) + C(x) * dC(x) / dxidxj )
    
    for(int j = 0 ; j < 3 ; j++)		// Kij's i, 'i' is index inside per triangle. Here j = 0 1 2 represents particle A B C
      for(int k = 0; k < 3; k++)		// Kij's j, 'j' is index inside per triangle. Here k = 0 1 2 represents particle A B C
      {
        for(int m = 0 ; m < 3; m++)			// loop inside per Kij 3*3 matrix
          for(int n = 0 ; n < 3; n++)		// loop inside per Kij 3*3 matrix
          {
            double identity = 0.0;
            double scale = 0.0;			// used as a scale factor for K-shear matrix
            if( m == n ) {identity = 1; scale = 1;}
            
            // K stretch
            double dFdz = area * dwudx[j] * dwudx[k] * (identity - wun[m] * wun[n]) * Cu / length_wu;				// (dCu2/dxm*dxn) * Cu
            dFdz += area * dwvdx[j] * dwvdx[k] * (identity - wvn[m] * wvn[n]) * Cv / length_wv;						// +(dCv2/dxm*dxn) * Cv ---> dC2/dxm*dxn * C
            dFdz += phiCstr[j][0][m] * phiCstr[k][0][n] + phiCstr[j][1][m] * phiCstr[k][1][n];						// +dC/dxm*transpose(dC/dxn)
            dFdz *= ( -groupTensileStiffness[group] );
            dFdz *= -1.0; // must flip because the internal force is on the left-hand side of equation
            
            // K shear
            double dFdz1 = area * Cs * scale * (dwudx[j] * dwvdx[k] + dwudx[k] * dwvdx[j]);
            dFdz1 += phiCsh[j][m] * phiCsh[k][n];
            //dFdz1 += phiCsh[j][m] * phiCsh[m][j];
            dFdz1 *= ( -groupShearStiffness[group] );
            dFdz1 *= -1.0; // must flip because the internal force is on the left-hand side of equation
            
            // Adding it to K topology
            int entry_row = 3 * triangles[3*i+j] + m; 
            int entry_col =	3 * inverseIndicesStretchAndShear[9*i+j*3+k] + n;
            K->AddEntry(entry_row, entry_col, dFdz);
            K->AddEntry(entry_row, entry_col, dFdz1);
          }
      }
  }
}

void ClothBW::AddBendStiffnessMatrix(double *u, SparseMatrix *K, int startQuad, int endQuad)
{
  for(int i = startQuad ; i < endQuad; i++)
  {
    // new ordering below
    int particleA = quadComponentIndices[6*i+0];	// index of non-edge vertex for triangle 1
    int particleB = quadComponentIndices[6*i+1];	// index of vertex 1 on the edge
    int particleC = quadComponentIndices[6*i+2];	// index of vertex 2 on the edge
    int particleD = quadComponentIndices[6*i+3];	// index of non-edge vertex for triangle 2
    
    int group1 = triangleGroups[quadComponentIndices[i*6+4]];			// group for first triangle
    int group2 = triangleGroups[quadComponentIndices[i*6+5]];			// group for second triangle
    
    double kBendU = (groupBendStiffnessU[group1] + groupBendStiffnessU[group2]) * 0.5;
    double kBendV = (groupBendStiffnessV[group1] + groupBendStiffnessV[group2]) * 0.5;
    double kBend = ( kBendU + kBendV ) * 0.5;
    
    // compute C(x)
    
    // first define some positions and accompanying vectors
    // triangle pictured below
    
    /* 
     1---------3
     / \       /
     /    \  5  /
     /  4   \   /
     /        \ /
     0----------2 
     */
    
    
    Vec3d positionA(&restPositions[3*particleA]);
    positionA += Vec3d(&u[3*particleA]);
    Vec3d positionB(&restPositions[3*particleB]);
    positionB += Vec3d(&u[3*particleB]);
    Vec3d positionC(&restPositions[3*particleC]);
    positionC += Vec3d(&u[3*particleC]);
    Vec3d positionD(&restPositions[3*particleD]);
    positionD += Vec3d(&u[3*particleD]);
    
    Vec3d vecBA = positionB - positionA;
    Vec3d vecCA = positionC - positionA;
    Vec3d vecBD = positionB - positionD;
    Vec3d vecCD = positionC - positionD;
    Vec3d vecBC = positionB - positionC;
    
    Vec3d vecCB = positionC - positionB;
    Vec3d vecAC = positionA - positionC;
    Vec3d vecDB = positionD - positionB;
    
    // now find normals for triangles 1 and 2
    Vec3d n1 = cross(vecCA, vecBA);
    Vec3d n2 = cross(vecBD, vecCD);
    Vec3d edge = vecBC;
    
    // now normalize the... normals
    Vec3d n1N = norm(n1);
    Vec3d n2N = norm(n2);
    Vec3d edgeNorm = norm(edge);
    
    // calculate sin theta and cos theta
    // sin theta = (norm1 x norm2) * edgeVector BC
    // cos theta = norm1 * norm2
    double sinTheta = dot(cross(n1N, n2N), edgeNorm);
    double cosTheta = dot(n1N, n2N);
    
    // calculate C(bend)
    double cbend = atan2(sinTheta, cosTheta);
    
    // account for resting angles?
    if (useRestAnglesForBendingForces)
      cbend -= restAngles[i];
    
//    // test this
//    while (cbend < -PI)
//      cbend += PI;
//    while (cbend > PI)
//      cbend -= PI;
//    // end test this
    
    double lengthN1 = len(n1);
    double lengthN2 = len(n2);
    double lengthE = len(edge);
    
    // calculate dn^a/dx, dn^b/dx, and dn^e/dx
    // dna/dx = Ss(q^a) from Pritchard's formula
    // dnb/dx = Ss(q^b)
    // dne/dx = I_s(q^e)
    
    // defined in Pritchard's paper
    // q^a = [BC][CA][BA][0]
    // q^b = [0] [CD][BD][BC]
    // q^e = [0] [1] [-1] [0]
    
    Vec3d qa[4];
    qa[0] = vecCB;
    qa[1] = vecAC;
    qa[2] = vecBA;
    qa[3][0] = qa[3][1] = qa[3][2] = 0.0;
    
    Vec3d qb[4];
    qb[0][0] = qb[0][1] = qb[0][2] =0.0;
    qb[1] = vecCD;
    qb[2] = vecDB;
    qb[3] = vecBC;
    
    /*
      dn1 is the first derivative of n1.

      n1 is a 3x1 vector, and d(n1)/d(vertexA) is a 3x3 matrix.
      Because we have 4 vertices, so dn1 is a 4x3x3 tensor.

      dn1[0] = d(n1)/d(vertexA) which is a 3x3 matrix
      dn1[1] = d(n1)/d(vertexB) which is a 3x3 matrix
      dn1[2] = d(n1)/d(vertexC) which is a 3x3 matrix
      dn1[3] = d(n1)/d(vertexD) which is a 3x3 matrix
      dn1[0][0] = d(n1)/d(vertexA[0]) which is a 3x1 vector
      (i.e., dn1[0][0] = [d(n1[0])/d(vertexA[0]), d(n1[1])/d(vertexA[0]), d(n1[2])/d(vertexA[0])] 
             dn1[0][1] = [d(n1[0])/d(vertexA[1]), d(n1[1])/d(vertexA[1]), d(n1[2])/d(vertexA[1])] 
             dn2[0][2] = [d(n1[0])/d(vertexA[2]), d(n1[1])/d(vertexA[2]), d(n1[2])/d(vertexA[2])] )

      dn1[derivative w.r.t. which vertex][d n1 component][d vertex component]
     */
    Vec3d dn1[4][3];
    Vec3d dn2[4][3];
    Vec3d dne[4][3];
    // each derivative conceptually structured as array of 3x3 matrices: 
    // [3x3 matrix] [3x3 matrix] [3x3 matrix] [3x3 matrix]
    for (int entry = 0; entry < 4; entry++)
    {
      //constructSkewSymmetric(qa[entry], dn1[entry], 1.0/lengthN1);
      //constructSkewSymmetric(qb[entry], dn2[entry], 1.0/lengthN2);
      
      double * dn1p = &(dn1[entry][0][0]);
      SKEW_MATRIX(qa[entry], dn1p);
      double invLengthN1 = 1.0 / lengthN1;
      for(int ii=0; ii<3; ii++)
        for(int jj=0; jj<3; jj++)
          dn1[entry][ii][jj] *= invLengthN1;
      
      double * dn2p = &(dn2[entry][0][0]);
      SKEW_MATRIX(qb[entry], dn2p);
      double invLengthN2 = 1.0 / lengthN2;
      for(int ii=0; ii<3; ii++)
        for(int jj=0; jj<3; jj++)
          dn2[entry][ii][jj] *= invLengthN2;
    }
    for (int outer=0; outer<3; outer++)
      for (int inner=0; inner<3; inner++)
      {
        dne[0][outer][inner] = 0.0;
        dne[3][outer][inner] = 0.0;
        
        if (outer == inner)
        {
          dne[1][outer][inner] = 1.0/lengthE;
          dne[2][outer][inner] = -1.0/lengthE;
        }
        else
        {
          dne[1][outer][inner] = 0.0;
          dne[2][outer][inner] = 0.0;
        }
      }

    
    /*
      sin and cos are scalars, so d(sin)/d(vertexA) is a 3x1 vector,
      and dSin is a 4x3x1 tensor

      dCos[0] = d(cosAngle)/d(vertexA) which is a 3x1 vector
      dCos[1] = d(cosAngle)/d(vertexB) which is a 3x1 vector
      dCos[2] = d(cosAngle)/d(vertexC) which is a 3x1 vector
      dCos[3] = d(cosAngle)/d(vertexD) which is a 3x1 vector

      dCos[derivative with respect to which vertex][and which x,y,z component]
     */
    // calculate dSin/dxm_s, dCos/dxm_s, and f = -k * C(x) * dC(x)/dx
    Vec3d dSin[4], dCos[4], dCb[4];

    for (int particle=0; particle<4; particle++)
      for (int xyz=0; xyz<3; xyz++)
      {
        dCos[particle][xyz] = dot(dn1[particle][xyz], n2N) + dot(n1N, dn2[particle][xyz]);
        
        dSin[particle][xyz] = dot ( cross(dn1[particle][xyz], n2N) + cross(n1N, dn2[particle][xyz]), edgeNorm );
        dSin[particle][xyz] += dot ( cross(n1N, n2N), dne[particle][xyz] );
        
        dCb[particle][xyz] = cosTheta * dSin[particle][xyz] - sinTheta * dCos[particle][xyz]; 
      }
    
    
    /*dqA=[0 -I  I  0;
     *     I  0  -I 0;
     *     -I I  0  0;
     *      0 0  0  0]
     
     *dqB=[0  0  0  0;
     *     0  0  I -I;
     *     0  -I 0  I;
     *     0  I  -I 0] */
    
    /*
      dQA[m][n][t][s], where (m,s) is a pair and (n,t) is a pair

      dQA = d(qA[m][s]) / d(x[n][t])

      (i.e., nt is the t-th component of the n-th vertex.
             For example, if n==1 and t==2, then it's the z-component of vertex 1.

             ms refers to the m-th row and s-th column of qA,
             like qA[m][s])

      dqA[0][0]/dx0x  dqA[0][0]/dx0y  dqA[0][0]/dx0z  dqA[0][0]/d1x ... dqA[0][0]/dx3z
      dqA[0][1]/dx0x  dqA[0][1]/dx0y  dqA[0][1]/dx0z  dqA[0][1]/d1x ... dqA[0][1]/dx3z
      ...
      dqA[3][2]/dx0x  dqA[3][2]/dx0y  dqA[3][2]/dx0z  dqA[3][2]/d1x ... dqA[3][2]/dx3z
     */
    double dQA[4][4][3][3], dQB[4][4][3][3];	// last 3 represents for the 3 value of a vector, dqA_m/dxn_t, dqB_m/dxn_t
    for(int m = 0; m < 4; m++)					//m loop
      for(int n = 0 ; n < 4; n++)				//n loop
        for(int x = 0; x < 3; x++)			//t loop
          for(int y = 0; y < 3; y++)		//3 values of a vector
          {
            dQA[m][n][x][y] = 0.0;
            dQB[m][n][x][y] = 0.0;
          }
    
    dQA[0][1][0][0] = -1.0; dQA[0][1][1][1] = -1.0; dQA[0][1][2][2] = -1.0;
    dQA[0][2][0][0] = 1.0;  dQA[0][2][1][1] = 1.0;  dQA[0][2][2][2] = 1.0;

    dQA[1][0][0][0] = 1.0;  dQA[1][0][1][1] = 1.0;  dQA[1][0][2][2] = 1.0;
    dQA[1][2][0][0] = -1.0; dQA[1][2][1][1] = -1.0; dQA[1][2][2][2] = -1.0;

    dQA[2][0][0][0] = -1.0; dQA[2][0][1][1] = -1.0; dQA[2][0][2][2] = -1.0;
    dQA[2][1][0][0] = 1.0;  dQA[2][1][1][1] = 1.0;  dQA[2][1][2][2] = 1.0;
    
    dQB[1][2][0][0] = 1.0;  dQB[1][2][1][1] = 1.0;  dQB[1][2][2][2] = 1.0;
    dQB[1][3][0][0] = -1.0; dQB[1][3][1][1] = -1.0; dQB[1][3][2][2] = -1.0;

    dQB[2][1][0][0] = -1.0; dQB[2][1][1][1] = -1.0; dQB[2][1][2][2] = -1.0;
    dQB[2][3][0][0] = 1.0;  dQB[2][3][1][1] = 1.0;  dQB[2][3][2][2] = 1.0;

    dQB[3][1][0][0] = 1.0;  dQB[3][1][1][1] = 1.0;  dQB[3][1][2][2] = 1.0;
    dQB[3][2][0][0] = -1.0; dQB[3][2][1][1] = -1.0; dQB[3][2][2][2] = -1.0;


    /*
      d2N1[m][n][s][t]
     */    
    Vec3d d2N1[4][4][3][3], d2N2[4][4][3][3];	// d2NA/dxm_s*dxn_t, d2N2/dxm_s*dxn_t
    
    //Vec3d d2N1Test[4][4][3][3], d2N2Test[4][4][3][3];
    
    for(int m = 0; m < 4; m++)					//m loop
      for(int n = 0 ; n < 4; n++)				//n loop
        for(int s = 0; s < 3; s++)			//s loop
          for(int t = 0; t < 3; t++)		//t loop
          {
            ClothBS_Sfunction(dQA[m][n][t], d2N1[m][n][s][t], s, 1.0/lengthN1);
            ClothBS_Sfunction(dQB[m][n][t], d2N2[m][n][s][t], s, 1.0/lengthN2);
          }
    
    /*
      d2Sin[m][n][s][t]
     */
    double d2Sin[4][4][3][3], d2Cos[4][4][3][3];
    for(int m = 0; m < 4; m++)					//m loop
      for(int n = 0 ; n < 4; n++)				//n loop
        for(int s = 0; s < 3; s++)			//s loop
          for(int t = 0; t < 3; t++)		//t loop
          {
            double temp_dCb = dot(d2N1[m][n][s][t], n2N); //changed to account for new indexing
            temp_dCb += dot(dn2[n][t], dn1[m][s]);
            temp_dCb += dot(dn1[n][t], dn2[m][s]);
            temp_dCb += dot(n1N, d2N2[m][n][s][t]); 

            ////////d2Cos[m][n][t][s] = temp_dCb; <- typo in Pritchard's note equation 31
            d2Cos[m][n][s][t] = temp_dCb; 
            
            Vec3d temp1, temp2;
            
            temp1 = cross (d2N1[m][n][s][t], n2N); 
            
            temp2 = temp1;
            
            temp1 = cross(dn1[m][s], dn2[n][t]);
            
            temp2 += temp1;
            
            temp1 = cross(dn1[n][t], dn2[m][s]);
            
            temp2 += temp1;
            
            temp1 = cross(n1N, d2N2[m][n][s][t]); 
            
            temp2 += temp1;	
            
            temp_dCb = dot(temp2, edgeNorm);
            
            temp1 = cross(dn1[m][s], n2N);
            
            temp2 = temp1;	
            
            temp1 = cross(n1N, dn2[m][s]);
            
            temp2 += temp1;	
            
            temp_dCb += dot(temp2, dne[n][t]);
            
            temp1 = cross(dn1[n][t], n2N);
            
            temp2 = temp1;
            
            temp1 = cross(n1N, dn2[n][t]);
            
            temp2 += temp1;
            
            temp_dCb += dot(temp2, dne[m][s]);
            
            d2Sin[m][n][s][t] = temp_dCb;
          }	
    
    for(int m = 0; m < 4; m++)					// Kij 'i
      for(int n = 0; n < 4; n++)				// Kij 'j
        for(int s = 0; s < 3; s++)			// 3*3 loop for K matrix block
          for(int t = 0; t < 3; t++)		// 3*3 loop for K matrix block
          {
            //double dFdz = -kBend * (dCb[m][s] * dCb[n][t] + cbend * (dCos[n][t] * dSin[m][s] + cosTheta * d2Sin[m][n][s][t] - dSin[n][t] * dCos[m][s] - sinTheta * d2Cos[m][n][s][t]) ); // changed to account for new indexing
            double dFdz = -kBend * ( dCb[m][s] * dCb[n][t] + cbend * (cosTheta * d2Sin[m][n][s][t] - sinTheta * d2Cos[m][n][s][t]) ); // observation that terms skew-symmetric terms are zero, Oct 2011
            dFdz *= -1.0; // must flip because the internal force is on the left-hand side of the equation
            int entry_row;
            int entry_col;
            // old order: CABD
            if     ( m == 0 )	entry_row = 3*particleA+s;
            else if( m == 1 )	entry_row = 3*particleB+s;
            else if( m == 2 )	entry_row = 3*particleC+s;
            else if( m == 3 )	entry_row = 3*particleD+s;
            // ------------------------------------------
            entry_col = 3 * inverseIndicesQuad[i*16+ m*4+ n] + t;
            K->AddEntry(entry_row, entry_col, dFdz);
          }   
  }
}

