/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "isotropic hyperelastic FEM" library , Copyright (C) 2014 USC         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin                            *
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

#include "isotropicHyperelasticFEM.h"
#include "matrixIO.h"
#include "mat3d.h"

#define SVD_singularValue_eps 1e-8

IsotropicHyperelasticFEM::IsotropicHyperelasticFEM(TetMesh * tetMesh_, IsotropicMaterial * isotropicMaterial_, double inversionThreshold_, bool addGravity_, double g_) :
  tetMesh(tetMesh_),
  isotropicMaterial(isotropicMaterial_),
  inversionThreshold(inversionThreshold_),
  addGravity(addGravity_), 
  g(g_)
{
  if (tetMesh->getNumElementVertices() != 4)
  {
    printf("Error in IsotropicHyperelasticFEM constructor: tets must have 4 vertices.\n");
    throw 1;
  }

  //printf("Entering IsotropicHyperelasticFEM::IsotropicHyperelasticFEM. Num vertices: %d. Num tets: %d.\n", tetMesh->getNumVertices(), tetMesh->getNumElements()); 

  int numElements = tetMesh->getNumElements();
  int numVertices = tetMesh->getNumVertices();

  // create space for F, U, Fhat, V, area weighted normals, and dmInverses (D_m^{-1})
  // each tet has numElementVertices vertices
  areaWeightedVertexNormals = (Vec3d*) malloc (sizeof(Vec3d) * numElements * tetMesh->getNumElementVertices());
  dmInverses = (Mat3d*) malloc (sizeof(Mat3d) * numElements);
  Fs = (Mat3d*) malloc (sizeof(Mat3d) * numElements);
  Fhats = (Vec3d*) malloc (sizeof(Vec3d)* numElements);
  Vs = (Mat3d*) malloc (sizeof(Mat3d) * numElements);
  Us = (Mat3d*) malloc (sizeof(Mat3d) * numElements);
  tetVolumes = (double*) malloc (sizeof(double) * numElements);

  currentVerticesPosition = (double*) malloc (sizeof(double) * 3 * numVertices);
  // save rest positions
  restVerticesPosition = (double*) malloc (sizeof(double) * 3 * numVertices);
  for (int i=0; i<numVertices; i++)
  {
    Vec3d * v = tetMesh->getVertex(i);
    restVerticesPosition[3*i+0] = (*v)[0];
    restVerticesPosition[3*i+1] = (*v)[1];
    restVerticesPosition[3*i+2] = (*v)[2];
  }

  ComputeTetVolumes();

  ComputeAreaWeightedVertexNormals(); //see p3 section 4 of [Irving 04]
  // precompute dmInverses (D_m^{-1}), which are needed to compute the 
  // deformation gradients at runtime ( F = D_s D_m^{-1} (see [Irving 04]) )
  PrepareDeformGrad(); //see p3 section 3 of [Irving 04]

  // build stiffness matrix skeleton
  // (i.e., create memory space for the non-zero entries of the stiffness matrix)
  SparseMatrix * stiffnessMatrixTopology;
  GetStiffnessMatrixTopology(&stiffnessMatrixTopology);

  // build acceleration indices so that we can quickly write the element stiffness matrices into the global stiffness matrix
  row_ = (int**) malloc (sizeof(int*) * numElements);
  column_ = (int**) malloc (sizeof(int*) * numElements);

  int numElementVertices = tetMesh->getNumElementVertices();

  for (int el=0; el < numElements; el++)
  {
    row_[el] = (int*) malloc (sizeof(int) * numElementVertices);
    column_[el] = (int*) malloc (sizeof(int) * numElementVertices * numElementVertices);

    for(int vertex=0; vertex<numElementVertices; vertex++)
      row_[el][vertex] = tetMesh->getVertexIndex(el, vertex);

    // seek for value row[j] in list associated with row[i]
    for(int i=0; i<numElementVertices; i++)
      for(int j=0; j<numElementVertices; j++)
        column_[el][numElementVertices * i + j] =
          stiffnessMatrixTopology->GetInverseIndex(3 * row_[el][i], 3 * row_[el][j]) / 3;
  }

  delete(stiffnessMatrixTopology);

  // DS = D_s
  // dDSdU = d D_s / d U
  // set dDSdU here, it's a constant matrix (does not change during the simulation)
  memset(dDSdU, 0.0, sizeof(double) * 108);
  dDSdU[tensor9x12Index(0,0,0,0)] = -1.0;
  dDSdU[tensor9x12Index(1,0,0,1)] = -1.0;
  dDSdU[tensor9x12Index(2,0,0,2)] = -1.0;
  dDSdU[tensor9x12Index(0,1,1,0)] = -1.0;
  dDSdU[tensor9x12Index(1,1,1,1)] = -1.0;
  dDSdU[tensor9x12Index(2,1,1,2)] = -1.0;
  dDSdU[tensor9x12Index(0,2,2,0)] = -1.0;
  dDSdU[tensor9x12Index(1,2,2,1)] = -1.0;
  dDSdU[tensor9x12Index(2,2,2,2)] = -1.0;
  dDSdU[tensor9x12Index(0,0,3,0)] = 1.0;
  dDSdU[tensor9x12Index(0,1,3,0)] = 1.0;
  dDSdU[tensor9x12Index(0,2,3,0)] = 1.0;
  dDSdU[tensor9x12Index(1,0,3,1)] = 1.0;
  dDSdU[tensor9x12Index(1,1,3,1)] = 1.0;
  dDSdU[tensor9x12Index(1,2,3,1)] = 1.0;
  dDSdU[tensor9x12Index(2,0,3,2)] = 1.0;
  dDSdU[tensor9x12Index(2,1,3,2)] = 1.0;
  dDSdU[tensor9x12Index(2,2,3,2)] = 1.0;

  dFdUs = (double*) malloc (sizeof(double) * 108 * numElements);
  Compute_dFdU(); // dF / dU is constant; precompute it

  // set the renumbering indices for conversion from Teran's order to row-major order
  rowMajorMatrixToTeran[0] = 0;
  rowMajorMatrixToTeran[1] = 3;
  rowMajorMatrixToTeran[2] = 5;
  rowMajorMatrixToTeran[3] = 4;
  rowMajorMatrixToTeran[4] = 1;
  rowMajorMatrixToTeran[5] = 7;
  rowMajorMatrixToTeran[6] = 6;
  rowMajorMatrixToTeran[7] = 8;
  rowMajorMatrixToTeran[8] = 2;

  for(int i=0; i<9; i++)
    teranToRowMajorMatrix[rowMajorMatrixToTeran[i]] = i;
}

IsotropicHyperelasticFEM::~IsotropicHyperelasticFEM()
{
  free(restVerticesPosition);
  free(currentVerticesPosition);
  free(Us);
  free(Vs);
  free(Fhats);
  free(Fs);
  free(dmInverses);
  free(areaWeightedVertexNormals);
  free(dFdUs);
  free(tetVolumes);

  int numElements = tetMesh->getNumElements();
  for (int el=0; el < numElements; el++)
  {
    free(row_[el]);
    free(column_[el]);
  }
  free(row_);
  free(column_);
}

/*
  Compute the elastic strain energy given the current vertex displacements u
*/
double IsotropicHyperelasticFEM::ComputeEnergy(double * u)
{
  int computationMode = COMPUTE_ENERGY;
  double energy;
  GetEnergyAndForceAndTangentStiffnessMatrixHelper(u, &energy, NULL, NULL, computationMode);
  return energy;
}

/*
  Compute the internal forces given the current vertex displacements u
*/
void IsotropicHyperelasticFEM::ComputeForces(double * u, double * internalForces)
{
  //printf("Entering IsotropicHyperelasticFEM::ComputeForces\n"); 
  int computationMode = IsotropicHyperelasticFEM::COMPUTE_INTERNALFORCES;
  GetEnergyAndForceAndTangentStiffnessMatrixHelper(u, NULL, internalForces, NULL, computationMode);
}

/*
  Create a matrix with the sparsity structure of the stiffness matrix.
  Note that this won't compute the actual values of the entries (they are set to zero).
  The sparsity structure does not change at runtime.
*/
void IsotropicHyperelasticFEM::GetStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
  int numVertices = tetMesh->getNumVertices();
  int numElementVertices = tetMesh->getNumElementVertices(); // This will always be 4.
  int * vertices = (int*) malloc (sizeof(int) * numElementVertices);

  // build the non-zero locations of the tangent stiffness matrix
  SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);
  int numElements = tetMesh->getNumElements();
  for (int el=0; el < numElements; el++)
  {
    for(int vertex=0; vertex<numElementVertices; vertex++)
      vertices[vertex] = tetMesh->getVertexIndex(el, vertex);

    for (int i=0; i<numElementVertices; i++)
      for (int j=0; j<numElementVertices; j++)
      {
        for(int k=0; k<3; k++)
          for(int l=0; l<3; l++)
          {
            // only add the entry if both vertices are free (non-fixed)
            // the corresponding elt is in row 3*i+k, column 3*j+l
            emptyMatrix->AddEntry( 3*vertices[i]+k, 3*vertices[j]+l, 0.0 );
          }
      }
  }

  *tangentStiffnessMatrix = new SparseMatrix(emptyMatrix);
  delete(emptyMatrix);

  free(vertices);
}

/*
  Get the tangent stiffness matrix given the current vertex displacement u
*/
void IsotropicHyperelasticFEM::GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix)
{
  int computationMode = COMPUTE_TANGENTSTIFFNESSMATRIX;
  GetEnergyAndForceAndTangentStiffnessMatrixHelper(u, NULL, NULL, tangentStiffnessMatrix, computationMode);
}

/*
  Get both internal forces and stiffness matrix, given the current vertex displacement u.
  Note: this is economic, because we have to compute the deformation gradients and do SVD on them only once.
*/
void IsotropicHyperelasticFEM::GetForceAndTangentStiffnessMatrix(double * u, double * internalForce, SparseMatrix * tangentStiffnessMatrix)
{
  int computationMode = COMPUTE_INTERNALFORCES | COMPUTE_TANGENTSTIFFNESSMATRIX;
  GetEnergyAndForceAndTangentStiffnessMatrixHelper(u, NULL, internalForce, tangentStiffnessMatrix, computationMode);
}

void IsotropicHyperelasticFEM::ComputeTetVolumes()
{
  int numElements = tetMesh->getNumElements();
  for (int el=0; el<numElements; el++)
    tetVolumes[el] = TetMesh::getTetVolume(tetMesh->getVertex(el, 0), tetMesh->getVertex(el, 1), tetMesh->getVertex(el, 2), tetMesh->getVertex(el, 3));
}

/*
  Compute the area-weighted vertex normals.
  See p3 section 4 of [Irving 04] for more details.
*/
void IsotropicHyperelasticFEM::ComputeAreaWeightedVertexNormals()
{
  int numElements = tetMesh->getNumElements();
  for (int el=0; el<numElements; el++)
  {
    Vec3d * va = tetMesh->getVertex(el, 0);
    Vec3d * vb = tetMesh->getVertex(el, 1);
    Vec3d * vc = tetMesh->getVertex(el, 2);
    Vec3d * vd = tetMesh->getVertex(el, 3);

    // compute normals for the four faces: acb, adc, abd, bcd
    Vec3d acbNormal = cross(*vc-*va, *vb-*va); 
    Vec3d adcNormal = cross(*vd-*va, *vc-*va); 
    Vec3d abdNormal = cross(*vb-*va, *vd-*va); 
    Vec3d bcdNormal = cross(*vc-*vb, *vd-*vb); 

    // if the tet vertices abcd form a positive orientation, the normals are now correct
    // otherwise, we need to flip them
    double orientation = dot(*vd-*va, cross(*vb-*va, *vc-*va));
    if (orientation < 0)
    {
      acbNormal *= -1.0;
      adcNormal *= -1.0;
      abdNormal *= -1.0;
      bcdNormal *= -1.0;
    }

    // triangle area = 0.5 |u x v|
    double acbArea = 0.5 * sqrt(dot(acbNormal, acbNormal));
    double adcArea = 0.5 * sqrt(dot(adcNormal, adcNormal));
    double abdArea = 0.5 * sqrt(dot(abdNormal, abdNormal));
    double bcdArea = 0.5 * sqrt(dot(bcdNormal, bcdNormal));

    // normalize
    acbNormal.normalize();
    adcNormal.normalize();
    abdNormal.normalize();
    bcdNormal.normalize();

    areaWeightedVertexNormals[4*el+0] = (acbArea * acbNormal + adcArea * adcNormal + abdArea * abdNormal) / 3.0;
    areaWeightedVertexNormals[4*el+1] = (acbArea * acbNormal + abdArea * abdNormal + bcdArea * bcdNormal) / 3.0;
    areaWeightedVertexNormals[4*el+2] = (acbArea * acbNormal + adcArea * adcNormal + bcdArea * bcdNormal) / 3.0;
    areaWeightedVertexNormals[4*el+3] = (adcArea * adcNormal + abdArea * abdNormal + bcdArea * bcdNormal) / 3.0;

    /*
      printf("--- areaWeightedVertexNormals ---\n");
      printf("a = "); areaWeightedVertexNormals[4*el+0].print();
      printf("b = "); areaWeightedVertexNormals[4*el+1].print();
      printf("c = "); areaWeightedVertexNormals[4*el+2].print();
      printf("d = "); areaWeightedVertexNormals[4*el+3].print();
    */
  }
}

/*
  Compute the inverse of the Dm matrices.
  The Dm is a 3x3 matrix where the columns are the edge vectors of a
  tet in rest configuration. See p3 section 3 of [Irving 04] for more details.
 */
void IsotropicHyperelasticFEM::PrepareDeformGrad()
{
  int numElements = tetMesh->getNumElements();
  for (int el=0; el<numElements; el++)
  {
    Vec3d * va = tetMesh->getVertex(el, 0);
    Vec3d * vb = tetMesh->getVertex(el, 1);
    Vec3d * vc = tetMesh->getVertex(el, 2);
    Vec3d * vd = tetMesh->getVertex(el, 3);

    Vec3d dm1 = *vd - *va;
    Vec3d dm2 = *vd - *vb;
    Vec3d dm3 = *vd - *vc;

    Mat3d tmp(dm1[0], dm2[0], dm3[0], dm1[1], dm2[1], dm3[1], dm1[2], dm2[2], dm3[2]);
    //printf("--- dm ---\n");
    //tmp.print();
    dmInverses[el] = inv(tmp);
    //printf("--- inv(dm) ---\n");
    //dmInverses[e].print();
  }
}

/*
  This function computes the energy, internal forces, and/or tangent stiffness matrix, as requested by the computationMode.
  It is declared virtual and is overloaded in the multicore (MT) derived class.
*/
int IsotropicHyperelasticFEM::GetEnergyAndForceAndTangentStiffnessMatrixHelper(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode)
{
  GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(u, energy, internalForces, tangentStiffnessMatrix, computationMode); // resets the energy, internal forces and/or tangent stiffness matrix to zero
  int code = GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(0, tetMesh->getNumElements(), u, energy, internalForces, tangentStiffnessMatrix, computationMode);
  return code;
}

// initializes the energy, internal forces, and/or stiffness matrix
void IsotropicHyperelasticFEM::GetEnergyAndForceAndTangentStiffnessMatrixHelperPrologue(double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode)
{
  // compute the current deformed positions
  int numVertices = tetMesh->getNumVertices();
  int numVertices3 = 3 * numVertices;
  for (int i=0; i<numVertices3; i++)
    currentVerticesPosition[i] = restVerticesPosition[i] + u[i];

  if (computationMode & COMPUTE_ENERGY)
    *energy = 0.0;

  if (computationMode & COMPUTE_INTERNALFORCES)
  {
    // reset internal forces
    if (addGravity)
    {
      for (int i=0; i<numVertices; i++)
      {
        internalForces[3*i+0] = 0.0;
        internalForces[3*i+1] = g; // gravity acts in negative-y direction; internal forces are opposite of external forces
        internalForces[3*i+2] = 0.0;
      }
    }
    else
    {
      // zero out the forces
      memset(internalForces, 0.0, sizeof(double) * numVertices3);
    }
  }

  if (computationMode & COMPUTE_TANGENTSTIFFNESSMATRIX)
  {
    // reset stiffness matrix
    tangentStiffnessMatrix->ResetToZero();
  }
}

/*
  This is the workhorse of the IFEM class. It computes strain energy,
  internal forces, and/or the tangent stiffness matrix for a subset of the elements, startEl <= el < endEl

  The strain energy is computed based on the user-specified material (e.g. StVK, neo-Hookean, Mooney-Rivlin)

  The internal forces computation is based on [Irving 04],
  and the stiffness matrix computation is based on 
  section 6 & 7 of [Teran 05].
*/
int IsotropicHyperelasticFEM::GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse(int startEl, int endEl, double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix, int computationMode)
{
  //printf("Entering IsotropicHyperelasticFEM::GetEnergyAndForceAndTangentStiffnessMatrixHelperWorkhorse\n"); 
  //printf("inversionThreshold=%G\n", inversionThreshold);

  int numElementVertices = tetMesh->getNumElementVertices();
  double energyResult = 0.0;
  //bool dropBelowThreshold = false; // becomes true when a principal stretch falls below the threshold; only used for printing out informative comments
  
  // traverse the elements and assemble strain energy, internal forces and tangent stiffness matrix
  int exitCode = 0;
  for (int el=startEl; el<endEl; el++)
  {
    /*
      Compute the deformation gradient F.
      F = Ds * inv(Dm), where Ds is a 3x3 matrix where
      the columns are edge vectors of a tet in the current deformation,
      and Dm is a 3x3 matrix where the columns are edge vectors of a tet in
      the rest configuration. See p3 section 3 of [Irving 04] for more details.
     */
    int vaIndex = 3 * tetMesh->getVertexIndex(el, 0);
    int vbIndex = 3 * tetMesh->getVertexIndex(el, 1);
    int vcIndex = 3 * tetMesh->getVertexIndex(el, 2);
    int vdIndex = 3 * tetMesh->getVertexIndex(el, 3);

    Vec3d va(currentVerticesPosition[vaIndex], currentVerticesPosition[vaIndex+1], currentVerticesPosition[vaIndex+2]);
    Vec3d vb(currentVerticesPosition[vbIndex], currentVerticesPosition[vbIndex+1], currentVerticesPosition[vbIndex+2]);
    Vec3d vc(currentVerticesPosition[vcIndex], currentVerticesPosition[vcIndex+1], currentVerticesPosition[vcIndex+2]);
    Vec3d vd(currentVerticesPosition[vdIndex], currentVerticesPosition[vdIndex+1], currentVerticesPosition[vdIndex+2]);

    Vec3d ds1 = vd - va;
    Vec3d ds2 = vd - vb;
    Vec3d ds3 = vd - vc;
    
    Mat3d tmp(ds1[0], ds2[0], ds3[0], ds1[1], ds2[1], ds3[1], ds1[2], ds2[2], ds3[2]);
    Fs[el] = tmp * dmInverses[el];
    //printf("F =\n");
    //Fs[el].print();

    /*
      The deformation gradient has now been computed and is available in Fs[el]
    */

    // perform modified SVD on the deformation gradient
    Mat3d & F = Fs[el];
    Mat3d & U = Us[el];
    Mat3d & V = Vs[el];
    Vec3d & Fhat = Fhats[el];
    int modifiedSVD = 1;
    if (SVD(F, U, Fhat, V, SVD_singularValue_eps, modifiedSVD) != 0)
    {
      printf("error in diagonalization, el=%d\n", el);
      exitCode = 1;
    }

    /*
      SVD for the deformation gradient has now been computed.
      It is available in Us[el], Fhats[el], Vs[el].
    */

    // clamp fHat if below the principal stretch threshold
    double fHat[3];
    int clamped = 0;
    for(int i = 0; i < 3; i++)
    {
      if(Fhats[el][i] < inversionThreshold)
      {
        //dropBelowThreshold = true;
        Fhats[el][i] = inversionThreshold;
        clamped |= (1 << i);
      }
    }
    fHat[0] = Fhats[el][0];
    fHat[1] = Fhats[el][1];
    fHat[2] = Fhats[el][2];
    clamped = 0; // disable clamping

    // query the user-provided isotropic material to compute the strain energy
    if (computationMode & COMPUTE_ENERGY)
      energyResult += tetVolumes[el] * ComputeEnergyFromStretches(el, fHat);

    if (computationMode & COMPUTE_INTERNALFORCES)
    {
      /*
        --- Now compute the internal forces ---
    
        The first Piola-Kirchhoff stress P is calculated by equation 1 
        in p3 section 5 of [Irving 04]. Once we have P, we can compute
        the nodal forces G=PBm as described in section 4 of [Irving 04]
      */

      double pHat[3];
      ComputeDiagonalPFromStretches(el, fHat, pHat); // calls the isotropic material to compute the diagonal P tensor, given the principal stretches in fHat
      Vec3d pHatv(pHat);
  
      // This is the 1st equation in p3 section 5 of [Irving 04]
      // P = Us[el] * diag(pHat) * trans(Vs[el])
      Mat3d P = Us[el];
      P.multiplyDiagRight(pHatv);
      P = P * trans(Vs[el]);

      //printf("--- P ---\n");
      //P.print();

      /*
        we compute the nodal forces by G=PBm as described in 
        section 4 of [Irving 04]
      */
      // multiply by 4 because each tet has 4 vertices
      Vec3d forceUpdateA = P * areaWeightedVertexNormals[4 * el + 0];
      Vec3d forceUpdateB = P * areaWeightedVertexNormals[4 * el + 1];
      Vec3d forceUpdateC = P * areaWeightedVertexNormals[4 * el + 2];
      Vec3d forceUpdateD = P * areaWeightedVertexNormals[4 * el + 3];
      // multiply by 3 because each force (at vertex) has 3 components
      int vIndexA = 3 * tetMesh->getVertexIndex(el,0);
      int vIndexB = 3 * tetMesh->getVertexIndex(el,1);
      int vIndexC = 3 * tetMesh->getVertexIndex(el,2);
      int vIndexD = 3 * tetMesh->getVertexIndex(el,3);

      internalForces[vIndexA+0] += forceUpdateA[0];
      internalForces[vIndexA+1] += forceUpdateA[1];
      internalForces[vIndexA+2] += forceUpdateA[2];
      
      internalForces[vIndexB+0] += forceUpdateB[0];
      internalForces[vIndexB+1] += forceUpdateB[1];
      internalForces[vIndexB+2] += forceUpdateB[2];
      
      internalForces[vIndexC+0] += forceUpdateC[0];
      internalForces[vIndexC+1] += forceUpdateC[1];
      internalForces[vIndexC+2] += forceUpdateC[2];
      
      internalForces[vIndexD+0] += forceUpdateD[0];
      internalForces[vIndexD+1] += forceUpdateD[1];
      internalForces[vIndexD+2] += forceUpdateD[2];
    }

    if (computationMode & COMPUTE_TANGENTSTIFFNESSMATRIX)
    {
      /*
        --- Now compute the tangent stiffness matrix ---

        This implementation is based on section 6 & 7 of [Teran 05].
        We go through each tet in the mesh and compute the element
        stiffness matrix K, and then we put the entries of K into
        the correct position of the final stiffness matrix (i.e.,
        the stiffness matrix for the entire mesh)
       */

      double K[144];
      ComputeTetK(el, K, clamped);

      // write matrices in place
      for(int vtxIndexA=0; vtxIndexA<4; vtxIndexA++)
        for(int vtxIndexB=0; vtxIndexB<4; vtxIndexB++)
        {
          int vtxA = tetMesh->getVertexIndex(el, vtxIndexA);
          //int vtxB = tetMesh->getVertexIndex(el, vtxIndexB);
          
          int columnIndexCompressed = column_[el][numElementVertices * vtxIndexA + vtxIndexB];
          
          for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
            {
              int row = 3 * vtxA + i;
              int columnIndex = 3 * columnIndexCompressed + j;
              double * value = &K[ELT(12, 3*vtxIndexA+i, 3*vtxIndexB+j)];
              
              tangentStiffnessMatrix->AddEntry(row, columnIndex, *value);
            }
        }
    }
  }

  //if (dropBelowThreshold)
  //  printf("Principal stretch dropped below %f\n", inversionThreshold);

  if (computationMode & COMPUTE_ENERGY)
    *energy = energyResult;

  return exitCode;
}

/*
  Converts a 3x3x3x4 tensor index to 9x12 matrix index

  i goes from [0, 2] inclusively
  j goes from [0, 2] inclusively
  m goes from [0, 3] inclusively
  n goes from [0, 2] inclusively
*/
int IsotropicHyperelasticFEM::tensor9x12Index(int i, int j, int m, int n)
{
  /*
    |  dF_00/du_v0x  dF_00/du_v0y  dF_00/du_v0z dF_00/du_v1x ...  dF_00/du_v3z  |
    |  dF_01/du_v0x  dF_01/du_v0y  dF_01/du_v0z dF_01/du_v1x ...  dF_01/du_v3z  |
    |  dF_02/du_v0x  dF_02/du_v0y  dF_02/du_v0z dF_02/du_v1x ...  dF_02/du_v3z  |
    |  dF_10/du_v0x  dF_10/du_v0y  dF_10/du_v0z dF_10/du_v1x ...  dF_10/du_v3z  |
    |                                      ...                                  |
    |  dF_22/du_v0x  dF_22/du_v0y  dF_22/du_v0z dF_22/du_v1x ...  dF_22/du_v3z  |    


                 | u_00 u_01 u_02 |   | v0x v0y v0z |
    where u is = | u_10 u_11 u_12 | = | v1x v1y v1z |
                 | u_20 u_21 u_22 |   | v2x v2y v2z |
                 | u_30 u_31 u_32 |   | v3x v3y v3z |

  */
  int rowIndex_in9x12Matrix = 3 * i + j;
  int columnIndex_in9x12Matrix = 3 * m + n;
  /*
    the resulting index is row major
    e.g.,
    | e[0]  e[1]  e[2] ...  e[9] |
    | e[10] e[11]      ...       |
  */
  return (12 * rowIndex_in9x12Matrix + columnIndex_in9x12Matrix);
}

/*
  Compute the derivative of the G (i.e., vertex forces) with respect to the
  deformation gradient F. The b0, b1, and b2 are the area weighted vertex 
  normal of the 1st, 2nd, and 3rd vertex of the tet. For example, b0=A0N0
  where A0 is the area and N is the vertex normal. For more information about
  the matrix G and b0,b1,b2, please see p3 section 4 of [Irving 04].

  Since G=PBm, so dG/dF = (dP/dF) * Bm
 */
void IsotropicHyperelasticFEM::Compute_dGdF(Vec3d * b0, Vec3d * b1, Vec3d * b2,
                                            double dPdF[81], double dGdF[81])
{
  //Both G and F are 3x3 matrices, so dGdF has 81 entries
    memset(dGdF, 0.0, sizeof(double) * 81);

  /*
           | ga_x gb_x gc_x |   | 0 1 2 |
    if G = | ga_y gb_y gc_y | = | 3 4 5 |
           | ga_z gb_z gc_z |   | 6 7 8 | 

    where ga, gb, gc are the nodal forces at vertex a,b,c

            | ba_0 bb_0 bc_0 |   | 0 1 2 |
    and B = | ba_1 bb_1 bc_1 | = | 3 4 5 |
            | ba_2 bb_2 bc_2 |   | 6 7 8 |

           | dga_x/dF_00 dga_x/dF_01 dga_x/dF_02 dga_x/dF_10 ... dga_x/dF_22 |
           | dga_y/dF_00 dga_y/dF_01 dga_y/dF_02 dga_y/dF_10 ... dga_y/dF_22 |
           | dga_z/dF_00 dga_z/dF_01 dga_z/dF_02 dga_z/dF_10 ... dga_z/dF_22 |
    dGdF = | dgb_x/dF_00 dgb_x/dF_01 dgb_x/dF_02 dgb_x/dF_10 ... dgb_x/dF_22 |
           |                                 ...                             |
           | dgc_z/dF_00 dgc_z/dF_01 dgc_z/dF_02 dgc_z/dF_10 ... dgc_z/dF_22 |
   */

  Vec3d * bVec[3] = { b0, b1, b2 };

  //dga_x/dF, dga_y/dF, dga_z/dF
  //dgb_x/dF, dgb_y/dF, dgb_z/dF
  //dgc_x/dF, dgc_y/dF, dgc_z/dF
  memset(dGdF, 0, sizeof(double) * 81);
  for(int abc=0; abc<3; abc++)
    for(int i=0; i<3; i++)
      for (int column=0; column<9; column++)
        for(int k=0; k<3; k++)
          dGdF[27 * abc + 9 * i + column] += dPdF[(3*i+k)*9+column]*(*(bVec[abc]))[k];

  /*
  printf("---- printing dGdF ----\n");
  for (int i=0; i<9; i++)
  {
    for (int j=0; j<9; j++)
      printf("%G ", dGdF[i*9+j]);
    printf("\n");
  }
  */
}

/*
  Compute the derivative of the deformation gradient F with respect
  to the displacement vector u.

  Because F = Ds*inv(Dm) (see p3 section 3 of [Irving 04]), we can
  compute dFdU as: dF/dU = (dDs/dU) * inv(Dm)

  The dF/dU is stored as dFdU in our code, and dDs/dU as dDSdU,
  and inv(Dm) as dmInv
 */
void IsotropicHyperelasticFEM::Compute_dFdU()
{
  int numElements = tetMesh->getNumElements();
  for (int el=0; el<numElements; el++)
  {
    double * dFdU = &dFdUs[108 * el];
    Mat3d & dmInv = dmInverses[el];
    for (int index=0; index<108; index++)
    {
      int n = index % 3;
      int m = (int)(index / 3) % 4;
      int j = (int)(index / 12) % 3;
      int i = (int)(index / 36) % 3;
      double result = 0.0;
      for (int k=0; k<3; k++)
	result += dDSdU[tensor9x12Index(i,k,m,n)] * dmInv[k][j];
      dFdU[tensor9x12Index(i,j,m,n)] = result;
    }
  }
}

/*
  Compute the element stiffness matrix K. 
  Since a tet has 12 dof (i.e., 4 vertices * 3 dof), the element K has 144 entries.

  The el is the index to the tet where we want to compute the K.

  We refer the reader to the first equation of p2 section 2 of the SCA 2011 poster
  "Invertible Isotropic Hyperelasticity using SVD Gradients"
  (the link can be found in Prof. Jernej Barbic's web page)
  for more information about how the element stiffness matrix K is computed
  from dP/dF (i.e., the derivative of the first Piola Kirchhoff stress P with
  respect to the deformation gradient F).

  The idea is that K = dG/du = (dG/dF)*(dF/du) = [(dP/dF)*Bm]*(dF/du).

  Parameter clamped tells the routine that the element is in the inversion handling
  regime where the singular values were clamped. The routine then computes correct
  K for this regime. I.e., if a singular value is below the clamping threshold, then
  altering this singular value has no effect on the internal forces, hence stiffness
  in that direction is zero.
  first bit of clamped: singular value 0 was clamped, 1=YES, 0=NO
  second bit of clamped: singular value 1 was clamped, 1=YES, 0=NO
  third bit of clamped: singular value 2 was clamped, 1=YES, 0=NO
 */
void IsotropicHyperelasticFEM::ComputeTetK(int el, double K[144], int clamped)
{
  /*
    dP/dF is a column major matrix, but is stored as a 1D vector
    
    | dP_11/dF_11  dP_11/dF_12  dP_11/dF_13  dP_11/dF_21 ... dP_11/dF_33 |
    | dP_12/dF_11  dP_12/dF_12  dP_12/dF_13  dP_12/dF_21 ... dP_12/dF_33 |
    |                              ...                                   |
    | dP_33/dF_11  dP_33/dF_12  dP_33/dF_13  dP_33/dF_21 ... dP_33/dF_33 |
  */
  double dPdF[81]; //in 9x9 matrix format
  double dGdF[81]; //in 9x9 matrix format

  Compute_dPdF(el, dPdF, clamped);
  Compute_dGdF(&(areaWeightedVertexNormals[4 * el + 0]), &(areaWeightedVertexNormals[4 * el + 1]),
               &(areaWeightedVertexNormals[4 * el + 2]), dPdF, dGdF);
  //dF_dU was already computed by the constructor before calling this function
  double * dFdU = &dFdUs[108 * el];

  // K is stored column-major (however, it doesn't matter because K is symmetric)
  for (int row=0; row<9; row++)
  {
    for (int column=0; column<12; column++)
    {
      double result = 0;
      for (int inner=0; inner<9; inner++)
      {
	//dGdF is 9x9, and dFdU is 9x12
	result += dGdF[9 * row + inner]*dFdU[12 * inner + column];
      }
      K[12 * column + row] = result;
    }
  }

  //The last three columns are combinations of the first nine columns.
  //The reason is that the nodal force of the 4th vertex equals to 
  //the minus of the sum of the 1st, 2nd, and 3rd vertices (see p3 
  //section 4 of [Irving 04]
  for (int row = 0; row < 12; row++)
  {
    //10th column
    K[12 * row +  9] = -K[12 * row + 0] - K[12 * row + 3] - K[12 * row + 6];
    //11th column
    K[12 * row + 10] = -K[12 * row + 1] - K[12 * row + 4] - K[12 * row + 7];
    //12th column
    K[12 * row + 11] = -K[12 * row + 2] - K[12 * row + 5] - K[12 * row + 8];
  }
}

double IsotropicHyperelasticFEM::ComputeEnergyFromStretches(int elementIndex, double * lambda)
{
  double invariants[3];

  double lambda2[3] = { lambda[0] * lambda[0], lambda[1] * lambda[1], lambda[2] * lambda[2] };
  double IC = lambda2[0] + lambda2[1] + lambda2[2];
  double IIC = lambda2[0] * lambda2[0] + lambda2[1] * lambda2[1] + lambda2[2] * lambda2[2];
  double IIIC = lambda2[0] * lambda2[1] * lambda2[2];

  invariants[0] = IC;
  invariants[1] = IIC;
  invariants[2] = IIIC;

  return isotropicMaterial->ComputeEnergy(elementIndex, invariants);
}

// compute diagonal Piola stress tensor from the three principal stretches
void IsotropicHyperelasticFEM::ComputeDiagonalPFromStretches(int elementIndex, double * lambda, double * PDiag)
{
  double invariants[3];

  double lambda2[3] = { lambda[0] * lambda[0], lambda[1] * lambda[1], lambda[2] * lambda[2] };
  double IC = lambda2[0] + lambda2[1] + lambda2[2];
  double IIC = lambda2[0] * lambda2[0] + lambda2[1] * lambda2[1] + lambda2[2] * lambda2[2];
  double IIIC = lambda2[0] * lambda2[1] * lambda2[2];

  invariants[0] = IC;
  invariants[1] = IIC;
  invariants[2] = IIIC;

  double dPsidI[3];
  
  isotropicMaterial->ComputeEnergyGradient(elementIndex, invariants, dPsidI);

  // PDiag = [ dI / dlambda ]^T * dPsidI

  double mat[9];
  mat[0] = 2.0 * lambda[0];
  mat[1] = 2.0 * lambda[1];
  mat[2] = 2.0 * lambda[2];
  mat[3] = 4.0 * lambda[0] * lambda[0] * lambda[0];
  mat[4] = 4.0 * lambda[1] * lambda[1] * lambda[1];
  mat[5] = 4.0 * lambda[2] * lambda[2] * lambda[2];
  mat[6] = 2.0 * lambda[0] * lambda2[1] * lambda2[2];
  mat[7] = 2.0 * lambda[1] * lambda2[0] * lambda2[2];
  mat[8] = 2.0 * lambda[2] * lambda2[0] * lambda2[1];

  Mat3d matM(mat);
  Vec3d dPsidIV(dPsidI);
  Vec3d result;

  result = trans(matM) * dPsidIV;
  result.convertToArray(PDiag);
}

/*
  The "i" goes from 0 to 2 inclusively
  The "j" goes from 0 to 2 inclusively
  See [Teran 05].
 */
inline double IsotropicHyperelasticFEM::gammaValue(int i, int j, double sigma[3], double invariants[3], double gradient[3], double hessian[6])
{
  /*
    The hessian is in order (11,12,13,22,23,33)
    | 11 12 13 |   | 0 1 2 |
    | 21 22 23 | = | 1 3 4 |
    | 31 32 33 |   | 2 4 5 |
  */

  double tempGammaVec1[3];
  tempGammaVec1[0] = 2.0 * sigma[i];
  tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
  tempGammaVec1[2] = 2.0 * invariants[2] / sigma[i];
  double tempGammaVec2[3];
  tempGammaVec2[0] = 2.0 * sigma[j];
  tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
  tempGammaVec2[2] = 2.0 * invariants[2] / sigma[j];
  double productResult[3];
  productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] + 
		      tempGammaVec2[2] * hessian[2]);
  productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] + 
		      tempGammaVec2[2] * hessian[4]);
  productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] + 
		      tempGammaVec2[2] * hessian[5]);
  return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
	  tempGammaVec1[2] * productResult[2] + 4.0 * invariants[2] * gradient[2] / (sigma[i] * sigma[j]));
}

// gradient of P with respect to F (9x9 matrix, row-major)
// see [Teran 05]
void IsotropicHyperelasticFEM::Compute_dPdF(int el, double dPdF[81], int clamped)
{
  double sigma[3] = { Fhats[el][0], Fhats[el][1], Fhats[el][2] };

  double sigma1square = sigma[0] * sigma[0];
  double sigma2square = sigma[1] * sigma[1];
  double sigma3square = sigma[2] * sigma[2];
  
  double invariants[3];
  invariants[0] = sigma1square + sigma2square + sigma3square;
  invariants[1] = (sigma1square * sigma1square + 
		   sigma2square * sigma2square +
		   sigma3square * sigma3square);
  invariants[2] = sigma1square * sigma2square * sigma3square;

  //double E[3];
  //E[0] = 0.5 * (Fhats[el][0] * Fhats[el][0] - 1);
  //E[1] = 0.5 * (Fhats[el][1] * Fhats[el][1] - 1);
  //E[2] = 0.5 * (Fhats[el][2] * Fhats[el][2] - 1);

  double gradient[3];
  isotropicMaterial->ComputeEnergyGradient(el, invariants, gradient);

  /*
    in order (11,12,13,22,23,33)
    | 11 12 13 |   | 0 1 2 |
    | 21 22 23 | = | 1 3 4 |
    | 31 32 33 |   | 2 4 5 |
  */
  double hessian[6];
  isotropicMaterial->ComputeEnergyHessian(el, invariants, hessian);

  // modify hessian to compute correct values if in the inversion handling regime
  if (clamped & 1) // first lambda was clamped (in inversion handling)
  {
    hessian[0] = hessian[1] = hessian[2] = 0.0;
  }

  if (clamped & 2) // second lambda was clamped (in inversion handling)
  {
    hessian[1] = hessian[3] = hessian[4] = 0.0;
  }

  if (clamped & 4) // third lambda was clamped (in inversion handling)
  {
    hessian[0] = hessian[1] = hessian[2] = hessian[4] = hessian[5] = 0.0;
  }

  double alpha11 = 2.0 * gradient[0] + 8.0 * sigma1square * gradient[1];
  double alpha22 = 2.0 * gradient[0] + 8.0 * sigma2square * gradient[1];
  double alpha33 = 2.0 * gradient[0] + 8.0 * sigma3square * gradient[1];
  double alpha12 = 2.0 * gradient[0] + 4.0 * (sigma1square+sigma2square) * gradient[1];
  double alpha13 = 2.0 * gradient[0] + 4.0 * (sigma1square+sigma3square) * gradient[1];
  double alpha23 = 2.0 * gradient[0] + 4.0 * (sigma2square+sigma3square) * gradient[1];

  double beta11 = 4.0 * sigma1square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma1square;
  double beta22 = 4.0 * sigma2square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma2square;
  double beta33 = 4.0 * sigma3square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma3square;
  double beta12 = 4.0 * sigma[0] * sigma[1] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[1]);
  double beta13 = 4.0 * sigma[0] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[2]);
  double beta23 = 4.0 * sigma[1] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[1] * sigma[2]);

  double gamma11 = gammaValue(0, 0, sigma, invariants, gradient, hessian);
  double gamma22 = gammaValue(1, 1, sigma, invariants, gradient, hessian);
  double gamma33 = gammaValue(2, 2, sigma, invariants, gradient, hessian);
  double gamma12 = gammaValue(0, 1, sigma, invariants, gradient, hessian);
  double gamma13 = gammaValue(0, 2, sigma, invariants, gradient, hessian);
  double gamma23 = gammaValue(1, 2, sigma, invariants, gradient, hessian);
  
  double x1111, x2222, x3333;
  double x2211, x3311, x3322;
  double x2121, x3131, x3232;
  double x2112, x3113, x3223;		 
   
  x1111 = alpha11 + beta11 + gamma11;
  x2222 = alpha22 + beta22 + gamma22;
  x3333 = alpha33 + beta33 + gamma33;

  x2211 = gamma12;
  x3311 = gamma13;
  x3322 = gamma23;

  x2121 = alpha12;
  x3131 = alpha13;
  x3232 = alpha23;

  x2112 = beta12;
  x3113 = beta13;
  x3223 = beta23;

  double dPdF_atFhat[81];
  memset(dPdF_atFhat, 0.0, sizeof(double) * 81);
  dPdF_atFhat[tensor9x9Index(0,0,0,0)] = x1111;
  dPdF_atFhat[tensor9x9Index(0,0,1,1)] = x2211;
  dPdF_atFhat[tensor9x9Index(0,0,2,2)] = x3311;

  dPdF_atFhat[tensor9x9Index(1,1,0,0)] = x2211;
  dPdF_atFhat[tensor9x9Index(1,1,1,1)] = x2222;
  dPdF_atFhat[tensor9x9Index(1,1,2,2)] = x3322;

  dPdF_atFhat[tensor9x9Index(2,2,0,0)] = x3311;
  dPdF_atFhat[tensor9x9Index(2,2,1,1)] = x3322;
  dPdF_atFhat[tensor9x9Index(2,2,2,2)] = x3333;

  dPdF_atFhat[tensor9x9Index(0,1,0,1)] = x2121;
  dPdF_atFhat[tensor9x9Index(0,1,1,0)] = x2112;

  dPdF_atFhat[tensor9x9Index(1,0,0,1)] = x2112;
  dPdF_atFhat[tensor9x9Index(1,0,1,0)] = x2121;

  dPdF_atFhat[tensor9x9Index(0,2,0,2)] = x3131;
  dPdF_atFhat[tensor9x9Index(0,2,2,0)] = x3113;

  dPdF_atFhat[tensor9x9Index(2,0,0,2)] = x3113;
  dPdF_atFhat[tensor9x9Index(2,0,2,0)] = x3131;

  dPdF_atFhat[tensor9x9Index(1,2,1,2)] = x3232;
  dPdF_atFhat[tensor9x9Index(1,2,2,1)] = x3223;

  dPdF_atFhat[tensor9x9Index(2,1,1,2)] = x3223;
  dPdF_atFhat[tensor9x9Index(2,1,2,1)] = x3232;

  /*
          | P_00 P_01 P_02 |        | F_00 F_01 F_02 |
    if P= | P_10 P_11 P_12 | and F= | F_10 F_11 F_12 |
          | P_20 P_21 P_22 |        | F_20 F_21 F_22 |

    | dP_00/dF_00  dP_00/dF_01 dP_00/dF_02 dP_00/dF_10 ... dP00/dF_22 |
    | dP_01/dF_00  dP_01/dF_01 dP_01/dF_02 dP_01/dF_10 ... dP01/dF_22 |
    | dP_02/dF_00  dP_02/dF_01 dP_02/dF_02 dP_02/dF_10 ... dP02/dF_22 |
    | dP_10/dF_00  dP_10/dF_01 dP_10/dF_02 dP_10/dF_10 ... dP10/dF_22 |
    |                               ...                               |
    | dP_22/dF_00  dP_22/dF_01 dP_22/dF_02 dP_22/dF_10 ... dP22/dF_22 |
   */

  Mat3d UT = trans(Us[el]); // trans(*U);
  Mat3d VT = trans(Vs[el]); // trans(*V);

  /*
    U->print();
    V->print();
    UT.print();
    VT.print();
  */

  double eiejVector[9];
  memset(eiejVector, 0.0, sizeof(double) * 9);
  memset(dPdF, 0.0, sizeof(double) * 81);
  for (int column=0; column<9; column++)
  {
    eiejVector[column] = 1.0;
    Mat3d ei_ej(eiejVector);
    Mat3d ut_eiej_v = UT*ei_ej*(Vs[el]);
    double ut_eiej_v_TeranVector[9]; //in Teran order
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v[0][0];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v[0][1];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v[0][2];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v[1][0];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[4]] = ut_eiej_v[1][1];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[5]] = ut_eiej_v[1][2];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[6]] = ut_eiej_v[2][0];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[7]] = ut_eiej_v[2][1];
    ut_eiej_v_TeranVector[rowMajorMatrixToTeran[8]] = ut_eiej_v[2][2];
    double dPdF_resultVector[9]; // not in Teran order
    for (int innerRow=0; innerRow<9; innerRow++)
    {
      double tempResult = 0.0;
      for (int innerColumn=0; innerColumn<9; innerColumn++)
      {
	tempResult += dPdF_atFhat[innerRow*9+innerColumn]*
	  ut_eiej_v_TeranVector[innerColumn];
      }
      dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
    }
    Mat3d dPdF_resultMatrix(dPdF_resultVector);
    Mat3d u_dpdf_vt = (Us[el])*dPdF_resultMatrix*VT;
    dPdF[column +  0] = u_dpdf_vt[0][0];
    dPdF[column +  9] = u_dpdf_vt[0][1];
    dPdF[column + 18] = u_dpdf_vt[0][2];
    dPdF[column + 27] = u_dpdf_vt[1][0];
    dPdF[column + 36] = u_dpdf_vt[1][1];
    dPdF[column + 45] = u_dpdf_vt[1][2];
    dPdF[column + 54] = u_dpdf_vt[2][0];
    dPdF[column + 63] = u_dpdf_vt[2][1];
    dPdF[column + 72] = u_dpdf_vt[2][2];
    // reset
    eiejVector[column] = 0.0;
  }

  /*
  printf("---- full dPdF ----\n");
  for (int i=0; i<9; i++)
  {
    for (int j=0; j<9; j++)
      printf("%G ", dPdF[i*9+j]);
    printf(";\n");
  }
  */
}

int IsotropicHyperelasticFEM::tensor9x9Index(int i, int j, int m, int n)
{
  /*
   |  dP_0/dF_0  dP_0/dF_4  dP_0/dF_8  ...  dP_0/dF_5  |
   |  dP_4/dF_0  dP_4/dF_4  dP_4/dF_8  ...  dP_4/dF_5  |
   |                         ...                       |
   |  dP_5/dF_0  dP_5/dF_4  dP_5/dF_8  ...  dP_5/dF_5  |
  */
  int rowIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * i + j];
  int columnIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * m + n];
  return (9 * rowIndex_in9x9Matrix + columnIndex_in9x9Matrix);
}

/*
  Compute damping forces based on the vertex velocities.
  See p6 section 6.2 in [Irving 04] for details.
*/
void IsotropicHyperelasticFEM::ComputeDampingForces(double dampingPsi, double dampingAlpha, double * u, double * uvel, double * dampingForces)
{
  Mat3d I(1.0); // identity matrix

  // --- damping forces ---
  int numElements = tetMesh->getNumElements();
  for (int el=0; el<numElements; el++)
  {
    // multiply by 3 because each velocity (at vertex) has 3 components
    int vaIndex = 3 * tetMesh->getVertexIndex(el, 0);
    int vbIndex = 3 * tetMesh->getVertexIndex(el, 1);
    int vcIndex = 3 * tetMesh->getVertexIndex(el, 2);
    int vdIndex = 3 * tetMesh->getVertexIndex(el, 3);

    Vec3d velocityA(uvel[vaIndex], uvel[vaIndex+1], uvel[vaIndex+2]);
    Vec3d velocityB(uvel[vbIndex], uvel[vbIndex+1], uvel[vbIndex+2]);
    Vec3d velocityC(uvel[vcIndex], uvel[vcIndex+1], uvel[vcIndex+2]);
    Vec3d velocityD(uvel[vdIndex], uvel[vdIndex+1], uvel[vdIndex+2]);

    Vec3d velocity1 = velocityD - velocityA;
    Vec3d velocity2 = velocityD - velocityB;
    Vec3d velocity3 = velocityD - velocityC;

    Mat3d tmp(velocity1[0], velocity2[0], velocity3[0], 
	      velocity1[1], velocity2[1], velocity3[1], 
	      velocity1[2], velocity2[2], velocity3[2]);

    Mat3d & U = Us[el];
    Mat3d & V = Vs[el];
    Mat3d FDotHat = trans(U) * (tmp * dmInverses[el]) * V;
    Mat3d Phat = 2 * dampingPsi * FDotHat + dampingAlpha * (FDotHat[0][0] + FDotHat[1][1] + FDotHat[2][2]) * I;
    Mat3d P = U * Phat * trans(V);

    Vec3d forceUpdate = P * areaWeightedVertexNormals[4 * el + 0];
    dampingForces[vaIndex+0] -= forceUpdate[0];
    dampingForces[vaIndex+1] -= forceUpdate[1];
    dampingForces[vaIndex+2] -= forceUpdate[2];

    forceUpdate = P * areaWeightedVertexNormals[4 * el + 1];
    dampingForces[vbIndex+0] -= forceUpdate[0];
    dampingForces[vbIndex+1] -= forceUpdate[1];
    dampingForces[vbIndex+2] -= forceUpdate[2];

    forceUpdate = P * areaWeightedVertexNormals[4 * el + 2];
    dampingForces[vcIndex+0] -= forceUpdate[0];
    dampingForces[vcIndex+1] -= forceUpdate[1];
    dampingForces[vcIndex+2] -= forceUpdate[2];

    forceUpdate = P * areaWeightedVertexNormals[4 * el + 3];
    dampingForces[vdIndex+0] -= forceUpdate[0];
    dampingForces[vdIndex+1] -= forceUpdate[1];
    dampingForces[vdIndex+2] -= forceUpdate[2];
  }
}

#undef modifiedSVD_singularValue_eps

