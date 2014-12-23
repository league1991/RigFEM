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

/*
  This class implements a mass-spring system in 3D.
  The mass spring system is given as a collection of particles (with masses),
  connected with springs (called "edges" in the code below).
  The springs can be organized into groups which share material 
  properties (stiffness and damping).
  The class can compute internal mass-spring system elastic forces,
  their gradient (the tangent stiffness matrix) and second derivative (Hessian).
  It can also compute the (diagonal) mass matrix of the structure.

  You can use the class in "massSpringSystemFromTetMesh.h" to create the mass-spring system directly from a tetrahedral mesh.
*/

#ifndef _MASS_SPRING_SYSTEM_H_
#define _MASS_SPRING_SYSTEM_H_

#include "sparseMatrix.h"

enum MassSpringSystemElementType {TET, CUBE};

class MassSpringSystem
{
public:

  // creates a mass spring system, by directly specifying the particle masses, and arbitrary springs connecting these masses
  // 'masses' is an array of length numParticles
  // 'restPositions' is an array of length 3xnumParticles
  // 'edges' is an integer array of length 2xnumEdges (giving integer indices of the two edge particles)
  // 'edgeGroups' is an integer array of length numEdges (giving the integer index of the material group to which each edge belongs)
  // 'groupStiffness' and 'groupDamping' are arrays that give the scalar stiffnesses and damping values for the edge groups
  // all indices in this class are 0-indexed
  MassSpringSystem(int numParticles, double * masses, double * restPositions,
    int numEdges, int * edges, int * edgeGroups, int numMaterialGroups, double * groupStiffness, double * groupDamping, int addGravity=0);

  // creates the mass spring system from a quad surface mesh (e.g., cloth)
  // each quad in the array 'quads' (total length 4xnumQuads) is given by 4 consecutive integer values
  // one group is created for the entire object
  // particle masses are computed automatically from the quad surface areas and the given surface mass density ('surfaceDensity' parameter)
  MassSpringSystem(int numParticles, double * restPositions, int numQuads, int * quads, double surfaceDensity, double tensileStiffness, double shearStiffness, double bendStiffness, double damping, int addGravity=0); 

  // creates the mass spring system from the elements of a 3D tetrahedral or cubic mesh
  // the tets or cubes are given in the array 'elements' (length 4x or 8x numElements) by 4 or 8 consecutive integer values
  // one group is created for the entire object
  // particle masses are computed automatically from the tet/cube volumes and the given mass density ('density' parameter): each tet gives 1/4 of its mass to each of its vertices, or each cube 1/8 of its mass
  MassSpringSystem(int numParticles, double * restPositions, MassSpringSystemElementType elementType, int numElements, int * elements, double density, double tensileStiffness, double damping, int addGravity=0); 

  // copy constructor
  MassSpringSystem(MassSpringSystem & massSpringSystem);

  virtual ~MassSpringSystem();

  void SetGravity(bool addGravity, double g=9.81) { this->addGravity = addGravity; this->g = g; } // if addGravity is enabled, ComputeForces will add the gravity force to the spring forces

  int GetNumParticles() { return numParticles; }
  int GetNumEdges() { return numEdges; }
  int * GetEdges() { return edges; }
  double * GetRestPositions() { return restPositions; }

  // creates the mass matrix (diagonal); each diagonal entry is expanded into a diagonal submatrix of size 'expanded' (typically, for 3D simulations, expanded should be 3)
  void GenerateMassMatrix(SparseMatrix ** M, int expanded=3);

  // compute the internal elastic force, under deformation u
  // important: the force f has the same sign as in the other deformable models in Vega, i.e., 
  //   it appears on the left side in equation M u'' + D u' + f = f_ext.
  //   If you want f to be interpreted as an external mass-spring force acting
  //   on the particles, you must flip the sign of f.
  //   Same comment applies to damping forces, stiffness matrices and their Hessian corrections.
  virtual void ComputeForce(double * u, double * f, bool addForce=false); // if addForce, f will be not be reset to zero prior to adding the forces

  // compute the damping force (it damps any relative velocities along each edge)
  virtual void ComputeDampingForce(double * uvel, double * f, bool addForce=false); 
  // compute the tangent stiffness matrix
  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); // call once to establish the location of sparse entries of the stiffness matrix
  virtual void ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix=false);
  // computes an approximation to dK, using the Hessian of internal forces, assuming the deformations change from u to u + du
  virtual void ComputeStiffnessMatrixCorrection(double * u, double * du, SparseMatrix * dK, bool addMatrix=false);

  // computes the gravitational force (result goes into f)
  void ComputeGravity(double * f, bool addForce=false);

  // creates a 3D triangle mesh, where each mass-spring system edge is one (degenerate) triangle
  // useful to visualize the mass-spring system
  // u is the deformation
  void CreateObjMesh(const char * filename, double * u = NULL); // if NULL, assumes zero deformation 

  // == advanced routines below ===

  void AddForce(double * u, double * f, int startEdge, int endEdge); 
  void AddStiffnessMatrix(double * u, SparseMatrix * K, int startEdge, int endEdge);
  void AddDampingForce(double * uvel, double * f, int startEdge, int endEdge); 
  void AddHessianApproximation(double * u, double * du, SparseMatrix * dK, int startEdge, int endEdge);

protected:

  friend class RenderSprings;

  void MassSpringSystemFromTets(int numParticles_, double * restPositions_, int numTets, int * tets, double density, double tensileStiffness, double damping);
  void MassSpringSystemFromCubes(int numParticles_, double * restPositions_, int numCubes, int * cubes, double density, double tensileStiffness, double damping);

  void GenerateMassSpringSystem(int numParticles, double * masses, double * restPositions, int numEdges, int * edges, int * edgeGroups, int numMaterialGroups, double * groupStiffness, double * groupDamping); // constructor helper function

  double GetTriangleSurfaceArea(double * p0, double * p1, double * p2);
  double GetTetVolume(double a[3], double b[3], double c[3], double d[3]);
  double GetCubeVolume(double a[3], double b[3], double c[3], double d[3], double e[3], double f[3], double g[3], double h[3]);

  int numParticles;
  double * masses;
  double * restPositions;
  double * restLengths;
  int numEdges;
  int * edges;
  int * inverseIndices;
  int * edgeGroups;

  int numMaterialGroups;
  double * groupStiffness, * groupDamping;

  int addGravity;
  double g;
};

#endif

