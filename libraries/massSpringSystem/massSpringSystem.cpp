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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <set>
#include "macros.h"
#include "massSpringSystem.h"
using namespace std;

MassSpringSystem::MassSpringSystem(int numParticles_, double * masses_, double * restPositions_, int numEdges_, int * edges_, int * edgeGroups_, int numMaterialGroups_, double * groupStiffness_, double * groupDamping_, int addGravity_) : addGravity(addGravity_), g(9.81)
{
  GenerateMassSpringSystem(numParticles_, masses_, restPositions_, numEdges_, edges_, edgeGroups_, numMaterialGroups_, groupStiffness_, groupDamping_);
}

MassSpringSystem::MassSpringSystem(int numParticles_, double * restPositions_, int numQuads, int * quads, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffness, double damping, int addGravity_): addGravity(addGravity_), g(9.81) // creates the mass spring system from a quad surface mesh (e.g., cloth)
{
  // set the masses
  double * masses_ = (double*) malloc (sizeof(double) * numParticles_);
  memset(masses_, 0, sizeof(double) * numParticles_); 
  for(int i=0; i<numQuads; i++)
  {
    // compute surface area of the quad
    int v0 = quads[4*i+0];
    int v1 = quads[4*i+1];
    int v2 = quads[4*i+2];
    int v3 = quads[4*i+3];
    double surfaceArea = GetTriangleSurfaceArea(&restPositions_[3*v0], &restPositions_[3*v1], &restPositions_[3*v2]) + GetTriangleSurfaceArea(&restPositions_[3*v2], &restPositions_[3*v0], &restPositions_[3*v3]);
    masses_[v0] += 0.25 * surfaceArea;
    masses_[v1] += 0.25 * surfaceArea;
    masses_[v2] += 0.25 * surfaceArea;
    masses_[v3] += 0.25 * surfaceArea;
  }

  // Generate springs:

  typedef pair<int,int> spring;
  #define SORTED(i,j) ( (i) <= (j) ? make_pair((i),(j)) : make_pair((j),(i)) )
  
  // first, all tensile springs (quad edges)
  set<spring> tensileSprings;
  for(int i=0; i<numQuads; i++)
  {
    int v0 = quads[4*i+0];
    int v1 = quads[4*i+1];
    int v2 = quads[4*i+2];
    int v3 = quads[4*i+3];
    tensileSprings.insert(SORTED(v0,v1));
    tensileSprings.insert(SORTED(v1,v2));
    tensileSprings.insert(SORTED(v2,v3));
    tensileSprings.insert(SORTED(v3,v0));
  }

  // then, all shear springs
  set<spring> shearSprings;
  for(int i=0; i<numQuads; i++)
  {
    int v0 = quads[4*i+0];
    int v1 = quads[4*i+1];
    int v2 = quads[4*i+2];
    int v3 = quads[4*i+3];
    shearSprings.insert(SORTED(v0,v2));
    shearSprings.insert(SORTED(v1,v3));
  }

  // and finally, all bend springs
  set<spring> bendSprings;
  // generate neighboring particles for every particle
  vector<set<int> > particleAxisNeighbors(numParticles_);
  vector<set<int> > particleNeighbors(numParticles_);
  for(int i=0; i<numQuads; i++)
  {
    int v0 = quads[4*i+0];
    int v1 = quads[4*i+1];
    int v2 = quads[4*i+2];
    int v3 = quads[4*i+3];

    // neighbors along each axis
    particleAxisNeighbors[v0].insert(v1);
    particleAxisNeighbors[v0].insert(v3);

    particleAxisNeighbors[v1].insert(v0);
    particleAxisNeighbors[v1].insert(v2);

    particleAxisNeighbors[v2].insert(v1);
    particleAxisNeighbors[v2].insert(v3);

    particleAxisNeighbors[v3].insert(v2);
    particleAxisNeighbors[v3].insert(v0);

    // particles sharing the same quad
    particleNeighbors[v0].insert(v0);
    particleNeighbors[v0].insert(v1);
    particleNeighbors[v0].insert(v2);
    particleNeighbors[v0].insert(v3);

    particleNeighbors[v1].insert(v0);
    particleNeighbors[v1].insert(v1);
    particleNeighbors[v1].insert(v2);
    particleNeighbors[v1].insert(v3);

    particleNeighbors[v2].insert(v0);
    particleNeighbors[v2].insert(v1);
    particleNeighbors[v2].insert(v2);
    particleNeighbors[v2].insert(v3);

    particleNeighbors[v3].insert(v0);
    particleNeighbors[v3].insert(v1);
    particleNeighbors[v3].insert(v2);
    particleNeighbors[v3].insert(v3);
  }
  // now, particleNeighbors[vtx] contains all particles that share a quad with vtx
  // particleAxisNeighbors[vtx] contains all particles that are neighbors of vtx in one of the axis directions (here, axes are aligned with the edge directions)
  // now, traverse all particles, and add bend neighbors:
  for(int i=0; i<numParticles_; i++)
  {
    for(set<int> :: iterator iter1 = particleAxisNeighbors[i].begin(); iter1 != particleAxisNeighbors[i].end(); iter1++)
      for(set<int> :: iterator iter2 = iter1; iter2 != particleAxisNeighbors[i].end(); iter2++)
      {
        if (particleNeighbors[*iter1].find(*iter2) == particleNeighbors[*iter1].end())
          bendSprings.insert(SORTED(*iter1, *iter2)); 
      }
  }
 
  int numEdges_ = tensileSprings.size() + shearSprings.size() + bendSprings.size();
  int * edges_ = (int*) malloc (sizeof(int) * 2 * numEdges_);
  int * edgeGroups_ = (int*) malloc (sizeof(int) * numEdges_);

  int count = 0;
  for(set<spring> :: iterator iter = tensileSprings.begin(); iter != tensileSprings.end(); iter++)
  {
    edges_[2*count+0] = iter->first;
    edges_[2*count+1] = iter->second;
    edgeGroups_[count] = 0; // tensile group
    count++;
  }

  for(set<spring> :: iterator iter = shearSprings.begin(); iter != shearSprings.end(); iter++)
  {
    edges_[2*count+0] = iter->first;
    edges_[2*count+1] = iter->second;
    edgeGroups_[count] = 1; // shear group
    count++;
  }

  for(set<spring> :: iterator iter = bendSprings.begin(); iter != bendSprings.end(); iter++)
  {
    edges_[2*count+0] = iter->first;
    edges_[2*count+1] = iter->second;
    edgeGroups_[count] = 2; // bend group
    count++;
  }

  printf("Num tensile springs: %d\n", (int)tensileSprings.size());
  printf("Num shear springs: %d\n", (int)shearSprings.size());
  printf("Num bend springs: %d\n", (int)bendSprings.size());

  int numMaterialGroups_ = 3;
  double * groupStiffness_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupStiffness_[0] = tensileStiffness;
  groupStiffness_[1] = shearStiffness;
  groupStiffness_[2] = bendStiffness;
  double * groupDamping_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupDamping_[0] = damping;
  groupDamping_[1] = damping;
  groupDamping_[2] = damping;

  GenerateMassSpringSystem(numParticles_, masses_, restPositions_, numEdges_, edges_, edgeGroups_, numMaterialGroups_, groupStiffness_, groupDamping_);

  free(masses_);
  free(edges_);
  free(edgeGroups_);
  free(groupStiffness_);
  free(groupDamping_);
}

MassSpringSystem::MassSpringSystem(int numParticles_, double * restPositions_, MassSpringSystemElementType elementType, int numElements, int * elements, double density, double tensileStiffness, double damping, int addGravity_) : addGravity(addGravity_)
{
  switch(elementType)
  {
    case TET:
      MassSpringSystemFromTets(numParticles_, restPositions_, numElements, elements, density, tensileStiffness, damping);
      break;
    case CUBE:
      MassSpringSystemFromCubes(numParticles_, restPositions_, numElements, elements, density, tensileStiffness, damping);
      break;
    default:
      printf("Unknown option\n");
      exit(1);
  }
}

// copy constructor
MassSpringSystem::MassSpringSystem(MassSpringSystem & massSpringSystem)
{
  numParticles = massSpringSystem.numParticles;

  masses = (double*) malloc (sizeof(double) * numParticles);
  memcpy(masses, massSpringSystem.masses, sizeof(double) * numParticles);

  restPositions = (double*) malloc (sizeof(double) * 3 * numParticles);
  memcpy(restPositions, massSpringSystem.restPositions, sizeof(double) * 3 * numParticles);

  numEdges = massSpringSystem.numEdges;

  restLengths = (double*) malloc (sizeof(double) * numEdges);
  memcpy(restLengths, massSpringSystem.restLengths, sizeof(double) * numEdges);

  edges = (int*) malloc (sizeof(int) * 2 * numEdges);
  memcpy(edges, massSpringSystem.edges, sizeof(int) * 2 * numEdges);

  inverseIndices = (int*) malloc (sizeof(int) * 4 * numEdges);
  memcpy(inverseIndices, massSpringSystem.inverseIndices, sizeof(int) * 4 * numEdges);

  edgeGroups = (int*) malloc (sizeof(int) * numEdges);
  memcpy(edgeGroups, massSpringSystem.edgeGroups, sizeof(int) * numEdges);

  numMaterialGroups = massSpringSystem.numMaterialGroups;

  groupStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupStiffness, massSpringSystem.groupStiffness, sizeof(double) * numMaterialGroups);

  groupDamping = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupDamping, massSpringSystem.groupDamping, sizeof(double) * numMaterialGroups);
  
  addGravity = massSpringSystem.addGravity;
  g = massSpringSystem.g;
}

MassSpringSystem::~MassSpringSystem()
{
  free(masses);
  free(restPositions);
  free(restLengths);
  free(edges);
  free(inverseIndices);
  free(edgeGroups);
  free(groupStiffness);
  free(groupDamping);
}

void MassSpringSystem::MassSpringSystemFromTets(int numParticles_, double * restPositions_, int numTets, int * tets, double density, double tensileStiffness, double damping)
{
  // set the masses
  double * masses_ = (double*) malloc (sizeof(double) * numParticles_);
  memset(masses_, 0, sizeof(double) * numParticles_); 
  for(int i=0; i<numTets; i++)
  {
    // compute surface area of the quad
    int v0 = tets[4*i+0];
    int v1 = tets[4*i+1];
    int v2 = tets[4*i+2];
    int v3 = tets[4*i+3];
    double tetVolume = GetTetVolume(&restPositions_[3*v0], &restPositions_[3*v1], &restPositions_[3*v2], &restPositions_[3*v3]);
    masses_[v0] += 0.25 * tetVolume;
    masses_[v1] += 0.25 * tetVolume;
    masses_[v2] += 0.25 * tetVolume;
    masses_[v3] += 0.25 * tetVolume;
  }

  // Generate springs:
  typedef pair<int,int> spring;
  #define SORTED(i,j) ( (i) <= (j) ? make_pair((i),(j)) : make_pair((j),(i)) )
  
  // first, all tensile springs (quad edges)
  set<spring> tensileSprings;
  for(int i=0; i<numTets; i++)
  {
    int v0 = tets[4*i+0];
    int v1 = tets[4*i+1];
    int v2 = tets[4*i+2];
    int v3 = tets[4*i+3];
    tensileSprings.insert(SORTED(v0,v1));
    tensileSprings.insert(SORTED(v0,v2));
    tensileSprings.insert(SORTED(v0,v3));

    tensileSprings.insert(SORTED(v1,v0));
    tensileSprings.insert(SORTED(v1,v2));
    tensileSprings.insert(SORTED(v1,v3));

    tensileSprings.insert(SORTED(v2,v0));
    tensileSprings.insert(SORTED(v2,v1));
    tensileSprings.insert(SORTED(v2,v3));

    tensileSprings.insert(SORTED(v3,v0));
    tensileSprings.insert(SORTED(v3,v1));
    tensileSprings.insert(SORTED(v3,v2));
  }

  int numEdges_ = tensileSprings.size(); 
  int * edges_ = (int*) malloc (sizeof(int) * 2 * numEdges_);
  int * edgeGroups_ = (int*) malloc (sizeof(int) * numEdges_);

  int count = 0;
  for(set<spring> :: iterator iter = tensileSprings.begin(); iter != tensileSprings.end(); iter++)
  {
    edges_[2*count+0] = iter->first;
    edges_[2*count+1] = iter->second;
    edgeGroups_[count] = 0; // tensile group
    count++;
  }

  printf("Num tensile springs: %d\n", (int)tensileSprings.size());

  int numMaterialGroups_ = 1;
  double * groupStiffness_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupStiffness_[0] = tensileStiffness;
  double * groupDamping_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupDamping_[0] = damping;

  GenerateMassSpringSystem(numParticles_, masses_, restPositions_, numEdges_, edges_, edgeGroups_, numMaterialGroups_, groupStiffness_, groupDamping_);

  free(masses_);
  free(edges_);
  free(edgeGroups_);
  free(groupStiffness_);
  free(groupDamping_);
}

void MassSpringSystem::MassSpringSystemFromCubes(int numParticles_, double * restPositions_, int numCubes, int * cubes, double density, double tensileStiffness, double damping)
{
  // set the masses
  double * masses_ = (double*) malloc (sizeof(double) * numParticles_);
  memset(masses_, 0, sizeof(double) * numParticles_); 
  for(int i=0; i<numCubes; i++)
  {
    // compute surface area of the quad
    int v0 = cubes[8*i+0];
    int v1 = cubes[8*i+1];
    int v2 = cubes[8*i+2];
    int v3 = cubes[8*i+3];
    int v4 = cubes[8*i+4];
    int v5 = cubes[8*i+5];
    int v6 = cubes[8*i+6];
    int v7 = cubes[8*i+7];
    double cubeVolume = GetCubeVolume(&restPositions_[3*v0], &restPositions_[3*v1], &restPositions_[3*v2], &restPositions_[3*v3],
                                      &restPositions_[3*v4], &restPositions_[3*v5], &restPositions_[3*v6], &restPositions_[3*v7]);
    masses_[v0] += 0.125 * cubeVolume;
    masses_[v1] += 0.125 * cubeVolume;
    masses_[v2] += 0.125 * cubeVolume;
    masses_[v3] += 0.125 * cubeVolume;
    masses_[v4] += 0.125 * cubeVolume;
    masses_[v5] += 0.125 * cubeVolume;
    masses_[v6] += 0.125 * cubeVolume;
    masses_[v7] += 0.125 * cubeVolume;
  }

  // Generate springs:
  typedef pair<int,int> spring;
  #define SORTED(i,j) ( (i) <= (j) ? make_pair((i),(j)) : make_pair((j),(i)) )
  
  // first, all tensile springs (quad edges)
  set<spring> tensileSprings;
  for(int i=0; i<numCubes; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      for(int k = j + 1; k < 8; k++)
      {
        tensileSprings.insert(SORTED(cubes[8*i+j], cubes[8*i+k]));
      }
    }
  }

  int numEdges_ = tensileSprings.size(); 
  int * edges_ = (int*) malloc (sizeof(int) * 2 * numEdges_);
  int * edgeGroups_ = (int*) malloc (sizeof(int) * numEdges_);

  int count = 0;
  for(set<spring> :: iterator iter = tensileSprings.begin(); iter != tensileSprings.end(); iter++)
  {
    edges_[2*count+0] = iter->first;
    edges_[2*count+1] = iter->second;
    edgeGroups_[count] = 0; // tensile group
    count++;
  }

  printf("Num tensile springs: %d\n", (int)tensileSprings.size());

  int numMaterialGroups_ = 1;
  double * groupStiffness_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupStiffness_[0] = tensileStiffness;
  double * groupDamping_ = (double*) malloc (sizeof(double) * numMaterialGroups_);
  groupDamping_[0] = damping;

  GenerateMassSpringSystem(numParticles_, masses_, restPositions_, numEdges_, edges_, edgeGroups_, numMaterialGroups_, groupStiffness_, groupDamping_);

  free(masses_);
  free(edges_);
  free(edgeGroups_);
  free(groupStiffness_);
  free(groupDamping_);
}

void MassSpringSystem::GenerateMassSpringSystem(int numParticles_, double * masses_, double * restPositions_, int numEdges_, int * edges_, int * edgeGroups_, int numMaterialGroups_, double * groupStiffness_, double * groupDamping_) 
{
  numParticles = numParticles_;
  numEdges = numEdges_;
  masses = (double*) malloc (sizeof(double) * numParticles);
  memcpy(masses, masses_, sizeof(double) * numParticles);
  restPositions = (double*) malloc (sizeof(double) * 3 * numParticles);
  memcpy(restPositions, restPositions_, sizeof(double) * 3 * numParticles);
  edges = (int*) malloc (sizeof(int) * 2 * numEdges);
  memcpy(edges, edges_, sizeof(int) * 2 * numEdges);

  edgeGroups = (int*) malloc (sizeof(int) * numEdges);
  memcpy(edgeGroups, edgeGroups_, sizeof(int) * numEdges);
  numMaterialGroups = numMaterialGroups_;
  groupStiffness = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupStiffness, groupStiffness_, sizeof(double) * numMaterialGroups);
  groupDamping = (double*) malloc (sizeof(double) * numMaterialGroups);
  memcpy(groupDamping, groupDamping_, sizeof(double) * numMaterialGroups);

  restLengths = (double*) malloc (sizeof(double) * numEdges);
  // compute rest lengths of springs
  for(int i=0; i<numEdges; i++)
  {
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    double restDisp[3];
    restDisp[0] = restPositions[3*particleB+0] - restPositions[3*particleA+0]; 
    restDisp[1] = restPositions[3*particleB+1] - restPositions[3*particleA+1];
    restDisp[2] = restPositions[3*particleB+2] - restPositions[3*particleA+2];

    restLengths[i] = sqrt(restDisp[0]*restDisp[0] + restDisp[1]*restDisp[1] + restDisp[2]*restDisp[2]);
  }

  // build inverse indices for stiffness matrix access
  SparseMatrixOutline skeletonOutline(numParticles);
  for(int i=0; i<numEdges; i++)
  {
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    skeletonOutline.AddEntry(particleA, particleA);
    skeletonOutline.AddEntry(particleA, particleB);
    skeletonOutline.AddEntry(particleB, particleA);
    skeletonOutline.AddEntry(particleB, particleB);
  }

  SparseMatrix skeleton(&skeletonOutline);
  inverseIndices = (int*) malloc (sizeof(int) * 4 * numEdges);
  for(int i=0; i<numEdges; i++)
  {
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];
    inverseIndices[4*i+0] = skeleton.GetInverseIndex(particleA, particleA);
    inverseIndices[4*i+1] = skeleton.GetInverseIndex(particleA, particleB);
    inverseIndices[4*i+2] = skeleton.GetInverseIndex(particleB, particleA);
    inverseIndices[4*i+3] = skeleton.GetInverseIndex(particleB, particleB);
  }
}

void MassSpringSystem::GenerateMassMatrix(SparseMatrix ** M, int expanded)
{
  SparseMatrixOutline outline(expanded * numParticles);
  for(int i=0; i<numParticles; i++)
    for(int j=0; j<expanded; j++)
      outline.AddEntry(expanded*i+j, expanded*i+j, masses[i]); 
  *M = new SparseMatrix(&outline);
}

double MassSpringSystem::GetTriangleSurfaceArea(double * p0, double * p1, double * p2)
{
  double s0[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
  double s1[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
  double crossp[3];
  CROSSPRODUCT(s0[0], s0[1], s0[2], s1[0], s1[1], s1[2], crossp[0], crossp[1], crossp[2]);
  return 0.5 * sqrt(crossp[0]*crossp[0] + crossp[1]*crossp[1] + crossp[2]*crossp[2]);
}

void MassSpringSystem::ComputeForce(double * u, double * f, bool addForce)
{
  if (!addForce)
    memset(f, 0, sizeof(double) * 3 * numParticles);

  AddForce(u, f, 0, numEdges);

  if (addGravity)
    ComputeGravity(f, true);
}

void MassSpringSystem::AddForce(double * u, double * f, int startEdge, int endEdge)
{
  for(int i=startEdge; i<endEdge; i++)
  {
    int group = edgeGroups[i];
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    double z[3]; // vector from A to B
    z[0] = restPositions[3*particleB+0] + u[3*particleB+0] - restPositions[3*particleA+0] - u[3*particleA+0];    
    z[1] = restPositions[3*particleB+1] + u[3*particleB+1] - restPositions[3*particleA+1] - u[3*particleA+1]; 
    z[2] = restPositions[3*particleB+2] + u[3*particleB+2] - restPositions[3*particleA+2] - u[3*particleA+2];    

    double len = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double force[3]; // force on particle A
    force[0] = groupStiffness[group] * (len - restLengths[i]) * z[0] / len;
    force[1] = groupStiffness[group] * (len - restLengths[i]) * z[1] / len;
    force[2] = groupStiffness[group] * (len - restLengths[i]) * z[2] / len;

    f[3*particleA+0] -= force[0];
    f[3*particleA+1] -= force[1];
    f[3*particleA+2] -= force[2];

    f[3*particleB+0] += force[0];
    f[3*particleB+1] += force[1];
    f[3*particleB+2] += force[2];
  }
}

void MassSpringSystem::GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology)
{
  SparseMatrixOutline KOutline(3*numParticles);

  for(int vtx=0; vtx<numParticles; vtx++)
  {
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        KOutline.AddEntry(3 * vtx + j, 3 * vtx + k);
  }

  for(int i=0; i<numEdges; i++)
  {
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

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

    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        KOutline.AddEntry(3 * particleA + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleA + j, 3 * particleB + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleA + k);
        KOutline.AddEntry(3 * particleB + j, 3 * particleB + k);
      }
  }

  *stiffnessMatrixTopology = new SparseMatrix(&KOutline);
}

void MassSpringSystem::ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix)
{
  if (!addMatrix)
    K->ResetToZero();

  AddStiffnessMatrix(u, K, 0, numEdges);
}

void MassSpringSystem::AddStiffnessMatrix(double * u, SparseMatrix * K, int startEdge, int endEdge)
{
  for(int i=startEdge; i<endEdge; i++)
  {
    int group = edgeGroups[i];
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    double dFdz[9];
    double z[3]; // z = rB - rA
    z[0] = restPositions[3*particleB+0] + u[3*particleB+0] - restPositions[3*particleA+0] - u[3*particleA+0];    
    z[1] = restPositions[3*particleB+1] + u[3*particleB+1] - restPositions[3*particleA+1] - u[3*particleA+1]; 
    z[2] = restPositions[3*particleB+2] + u[3*particleB+2] - restPositions[3*particleA+2] - u[3*particleA+2];    
     
    double len = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double invLen = 1.0 / len;
    z[0] *= invLen;
    z[1] *= invLen;
    z[2] *= invLen;
    memset(dFdz, 0, sizeof(double) * 9);
    dFdz[0] = 1.0 - restLengths[i] * invLen;
    dFdz[4] = 1.0 - restLengths[i] * invLen;
    dFdz[8] = 1.0 - restLengths[i] * invLen;
    
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        dFdz[3*k+j] += restLengths[i] * z[j] * z[k] * invLen;

    for(int j=0; j<9; j++)
      dFdz[j] *= groupStiffness[group];

    // write matrices in place
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        K->AddEntry(3 * particleA + j, 3 * inverseIndices[4*i+0] + k, +dFdz[3*k+j]);
        K->AddEntry(3 * particleA + j, 3 * inverseIndices[4*i+1] + k, -dFdz[3*k+j]);
        K->AddEntry(3 * particleB + j, 3 * inverseIndices[4*i+2] + k, -dFdz[3*k+j]);
        K->AddEntry(3 * particleB + j, 3 * inverseIndices[4*i+3] + k, +dFdz[3*k+j]);
      }
  }
}

void MassSpringSystem::ComputeDampingForce(double * uvel, double * f, bool addForce)
{
  if (!addForce)
    memset(f, 0, sizeof(double) * 3 * numParticles);

  AddDampingForce(uvel, f, 0, numEdges);
}

void MassSpringSystem::AddDampingForce(double * uvel, double * f, int startEdge, int endEdge)
{
  for(int i=startEdge; i<endEdge; i++)
  {
    int group = edgeGroups[i];
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    double z[3]; // z = rB' - rA';
    z[0] = uvel[3*particleB+0] - uvel[3*particleA+0];    
    z[1] = uvel[3*particleB+1] - uvel[3*particleA+1];    
    z[2] = uvel[3*particleB+2] - uvel[3*particleA+2];    
     
    double force[3]; // damping force on particle A
    force[0] = groupDamping[group] * z[0];
    force[1] = groupDamping[group] * z[1];
    force[2] = groupDamping[group] * z[2];

    f[3*particleA+0] -= force[0];
    f[3*particleA+1] -= force[1];
    f[3*particleA+2] -= force[2];

    f[3*particleB+0] += force[0];
    f[3*particleB+1] += force[1];
    f[3*particleB+2] += force[2];
  }
}

void MassSpringSystem::ComputeStiffnessMatrixCorrection(double * u, double * du, SparseMatrix * dK, bool addMatrix)
{
  if (!addMatrix)
    dK->ResetToZero();
  AddHessianApproximation(u, du, dK, 0, numEdges);
}

void MassSpringSystem::AddHessianApproximation(double * u, double * du, SparseMatrix * dK, int startEdge, int endEdge)
{
  for(int i=startEdge; i<endEdge; i++)
  {
    int group = edgeGroups[i];
    int particleA = edges[2*i+0];
    int particleB = edges[2*i+1];

    double z[3]; // z = rB - rA
    z[0] = restPositions[3*particleB+0] + u[3*particleB+0] - restPositions[3*particleA+0] - u[3*particleA+0];    
    z[1] = restPositions[3*particleB+1] + u[3*particleB+1] - restPositions[3*particleA+1] - u[3*particleA+1]; 
    z[2] = restPositions[3*particleB+2] + u[3*particleB+2] - restPositions[3*particleA+2] - u[3*particleA+2];    
    double len = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);

    double invlen = 1.0 / len;
    double invlen3 = invlen * invlen * invlen;
    double core[9];
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        core[k*3+j] = - 3.0 * invlen3 * z[j] * z[k];
    core[0] += 1.0 * invlen;
    core[4] += 1.0 * invlen;
    core[8] += 1.0 * invlen;

    double xA[3];
    xA[0] = du[3*particleA+0];
    xA[1] = du[3*particleA+1];
    xA[2] = du[3*particleA+2];
    double dAdz_xA[9];
    memset(dAdz_xA, 0, sizeof(double) * 9);
    double factorA = (z[0]*xA[0] + z[1]*xA[1] + z[2]*xA[2]) * invlen * invlen;
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        dAdz_xA[k*3+j] = factorA * core[k*3+j] + invlen3 * z[j]*xA[k] + invlen3 * xA[j]*z[k];
    
    double xB[3];
    xB[0] = du[3*particleB+0];
    xB[1] = du[3*particleB+1];
    xB[2] = du[3*particleB+2];
    double dAdz_xB[9];
    memset(dAdz_xB, 0, sizeof(double) * 9);
    double factorB = (z[0]*xB[0] + z[1]*xB[1] + z[2]*xB[2])/(len*len);
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        dAdz_xB[k*3+j] = factorB * core[k*3+j] + invlen3 * z[j]*xB[k] + invlen3 * xB[j]*z[k];
    
    for(int j=0; j<9; j++)
    {
      dAdz_xA[j] *= groupStiffness[group] * restLengths[i];
      dAdz_xB[j] *= groupStiffness[group] * restLengths[i];
    }

    double block[9];
    for(int j=0; j<9; j++)
      block[j] = dAdz_xA[j] - dAdz_xB[j];

    // write matrices in place
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        dK->AddEntry(3 * particleA + j, 3 * inverseIndices[4*i+0] + k, -block[3*k+j]);
        dK->AddEntry(3 * particleA + j, 3 * inverseIndices[4*i+1] + k, +block[3*k+j]);
        dK->AddEntry(3 * particleB + j, 3 * inverseIndices[4*i+2] + k, +block[3*k+j]);
        dK->AddEntry(3 * particleB + j, 3 * inverseIndices[4*i+3] + k, -block[3*k+j]);
      }
  }
}

void MassSpringSystem::ComputeGravity(double * f, bool addGravity)
{
  if (!addGravity)
    memset(f, 0, sizeof(double) * 3 * numParticles);

  for(int i=0; i<numParticles; i++)
    f[3*i+1] += -g * masses[i];
}

double MassSpringSystem::GetTetVolume(double a[3], double b[3], double c[3], double d[3])
{
  // volume = 1/6 * | (a-d) . ((b-d) x (c-d)) |

  double x[3] = { a[0] - d[0], a[1] - d[1], a[2] - d[2] };
  double y[3] = { b[0] - d[0], b[1] - d[1], b[2] - d[2] };
  double z[3] = { c[0] - d[0], c[1] - d[1], c[2] - d[2] };

  double w[3];
  CROSSPRODUCT(y[0], y[1], y[2], z[0], z[1], z[2], w[0], w[1], w[2]);
  
  return (1.0 / 6 * fabs(x[0] * w[0] + x[1] * w[1] + x[2] * w[2]));
}

double MassSpringSystem::GetCubeVolume(double a[3], double b[3], double c[3], double d[3], double e[3], double f[3], double g[3], double h[3])
{
  //all 8 verts are passed in case the volume method is changed later
  return (g[0] - a[0]) * (g[1] - a[1]) * (g[2] - a[2]);
}

void MassSpringSystem::CreateObjMesh(const char * filename, double * u)
{
  FILE * fout = fopen(filename, "w");

  if (u != NULL)
    for(int i=0; i<numParticles; i++)
      fprintf(fout, "v %G %G %G\n", 
        restPositions[3*i+0] + u[3*i+0], 
        restPositions[3*i+1] + u[3*i+1],
        restPositions[3*i+2] + u[3*i+2]);
  else
    for(int i=0; i<numParticles; i++)
      fprintf(fout, "v %G %G %G\n", 
        restPositions[3*i+0], restPositions[3*i+1], restPositions[3*i+2]);

  fprintf(fout, "\n");
  fprintf(fout, "vn 1 0 0\n");
  fprintf(fout, "\n");

  for(int i=0; i<numEdges; i++)
    fprintf(fout, "f %d//1 %d//1 %d//1\n", edges[2*i+0] + 1, edges[2*i+0] + 1, edges[2*i+1] + 1);

  fclose(fout);
}

