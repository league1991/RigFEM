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

/*
 This class implements the cloth model from [Baraff and Witkin 1998].
 It can compute both the internal elastic forces their gradient (tangent stiffness matrix).
 It computes stretch, shear, and bend forces.
 Damping is not implemented, but you can use the damping in the Vega integrator class.
*/

#ifndef _CLOTHBW_H_
#define _CLOTHBW_H_

#include "sparseMatrix.h"

class ClothBW
{
public:
  // creates the cloth elastic model, from a given triangle mesh
  // "masses" is an array of length "numParticles"
  // "restPositions" is an array of length 3x"numParticles"
  // "triangles" is an integer array of length 3x"numTriangles" (giving integer indices of the three particle forming a triangle)
  // "triangleGroups" is an integer array of length "numTriangles" (giving the integer index of the material group to which each triangle belongs)
  // "groupTensileStiffness", "groupShearStiffness", "groupBendStiffness", and "groupDamping" are arrays that give the scalar stiffness and damping parameters for each material groups
  // "triangleUVs" is an double array of length 3x2x"numTriangles", indicating the uv for every vertex; this array can be user-provided, or the constructor can compute it automatically
  // all indices in this class are 0-indexed
  // constructor that does not require triangleUVs input (it computes UVs automatically; note: the UVs are continuous only within each triangle; the UV map is not global (which is fine, provided one does not want to simulate anisotropic effects) )
  ClothBW(int numParticles, double * masses, double * restPositions, int numTriangles, int * triangles, int * triangleGroups, int numMaterialGroups, double * groupTensileStiffness, double * groupShearStiffness, double * groupBendStiffnessU, double * groupBendStiffnessV, double * groupDamping, int addGravity=0);
  
  // constructor with triangleUVs input
  ClothBW(int numParticles, double * masses, double * restPositions, int numTriangles, int * triangles, double * triangleUVs, int * triangleGroups, int numMaterialGroups, double * groupTensileStiffness, double * groupShearStiffness, double * groupBendStiffnessU, double * groupBendStiffnessV, double * groupDamping, int addGravity=0);	
    
  // copy constructor
  ClothBW(ClothBW & ClothBW);
  virtual ~ClothBW();
    
  void SetRestUVStretchValues(double bu, double bv) { this->bu = bu; this->bv = bv; } // these are the b_u and b_v stretch values from Equation (10) in the BW paper (default values are 1.0)
    
  int GetNumParticles() { return numParticles; }
  int GetNumTriangles() { return numTriangles; }
  int * GetTriangles() { return triangles; }
  double * GetRestPositions() { return restPositions;  }
    
  // creates the mass matrix (which is diagonal); each diagonal entry is expanded into a diagonal submatrix of size 'expanded' (typically, for 3D simulations, expanded should be 3)
  void GenerateMassMatrix(SparseMatrix ** M, int expanded=3);

  // computes the gravitational force (result goes into f)
  void SetGravity(bool addGravity, double g=9.81) { this->addGravity = addGravity; this->g = g; } // if addGravity is enabled, ComputeForces will add the gravity force to the elastic forces
  void ComputeGravity(double * f, bool addForce=false);
    
  // === compute elastic forces ===    
    
  // compute the internal elastic force, under deformation u
  // note: the force has the sign of the left side of the dynamic equation, Mu'' + Du' + f_int(u) = f_ext(t), i.e., f_int(u), that is, **opposite** to an external force f_ext(t) acting on the body 
  virtual void ComputeForce(double * u, double * f, bool addForce=false); // if addForce is "true", f will be not be reset to zero prior to adding the forces
  // compute the tangent stiffness matrix of the elastic force
  void GenerateStiffnessMatrixTopology(SparseMatrix ** K); // call once to establish the location of sparse entries of the stiffness matrix
  virtual void ComputeStiffnessMatrix(double * u, SparseMatrix * K, bool addMatrix=false);
    
  // compute the damping force
  // unimplemented
  // (use damping provided in the integrator class)
  virtual void ComputeDampingForce(double * u, double * uvel, double * f, bool addForce=false); 

  // this allows users to toggle on/off computation of the stretch/shear forces, bend forces,
  // stretch/shear stiffness matrix, and bend stiffness matrix
  // mode[0] = computeStretchAndShearForce
  // mode[1] = computeBendForce
  // mode[2] = computeStretchAndShearStiffnessMatrices
  // mode[3] = computeBendStiffnessMatrices
  void SetComputationMode(bool mode[4]);
    
  // allows user to toggle on/off the use of rest angles in bend force/stiffness matrix
  // calculations; if set to 1, bend force/matrix will be computed in relation to the quad's rest
  // angle. if set to 0, bend force/matrix will be computed in relation to a flat angle of 0.
  // default is 1.
  void UseRestAnglesForBendingForces(int useRestAnglesForBend)
  { 
    this->useRestAnglesForBendingForces = useRestAnglesForBend;
  }

  // === routines that only compute one type of a force ===
    
  void AddForce(const double * u, double * f, int startTriangle, int endTriangle, int startQuad, int endQuad);
  void AddStretchAndShearForce(const double * u, double * f, int startTriangle, int endTriangle);
  void AddBendForce(const double * u, double * f, int startQuad, int endQuad);
    
  void AddStiffnessMatrix(double * u, SparseMatrix * K, int startTriangle, int endTriangle, int startQuad, int endQuad);
    
  void AddStretchAndShearStiffnessMatrix(double * u, SparseMatrix * K, int startTriangle, int endTriangle);
  void AddBendStiffnessMatrix(double * u, SparseMatrix * K, int startQuad, int endQuad);
    
  // unimplemented
  // (use damping provided in the integrator class)
  void AddDampingForce(double * uvel, double * f, int startTriangle, int endTriangle);
        
protected:
  double GetTriangleSurfaceArea(double * p0, double * p1, double * p2);
  void GenerateBW(int numParticles, double * masses, double * restPositions, int numTriangles, int * triangles, double * triangleUVs, int * triangleGroups, int numMaterialGroups, double * groupTensileStiffness, double * groupShearStiffness, double * groupBendStiffnessU, double * groupBendStiffnessV, double * groupDamping, int addGravity=0);
    
  int numParticles;
  double * masses;
  double * restPositions;
  int numTriangles;
  int * triangles;
  double * triangleUVs;
  int * inverseIndicesStretchAndShear;
  int * triangleGroups;
    
  // Internal variables for bending computation:
  // "numQuads" is number of 'quads' made up of two adjacent triangles. Important because each quad contains a bendable edge (the edge shared by the two triangles).
  // "restAngles" is array of size numQuads containing resting angles of edges that can bend (indexed by numQuads)
  // "inverseIndicesQuad" is numQuads x 16 interger array
  // "quadComponentIndices" is numQuads x 6 interger array
  int numQuads;
  double * restAngles;
  int * inverseIndicesQuad;
  int * quadComponentIndices;
    
  int numMaterialGroups;
  double * groupTensileStiffness; 
  double * groupShearStiffness;
  double * groupBendStiffnessU;
  double * groupBendStiffnessV;
  double * groupDamping;
    
  double bu; // stretch constraint u (default = 1.0)
  double bv; // stretch constraint v (default = 1.0)
    
  int addGravity;
  double g;
    
  bool _computationConditions[4];
  int useRestAnglesForBendingForces;
};

#endif

