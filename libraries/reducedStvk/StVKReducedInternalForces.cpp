/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.0                               *
 *                                                                       *
 * "reducedStvk" library , Copyright (C) 2007 CMU, 2009 MIT              *
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

#include "lapack-headers.h"
#include "matrixIO.h"
#if defined(WIN32) || defined(linux)
  #include "mkl_service.h"
#endif
#include "matrixMacros.h"
#include "matrixProjection.h"
#include "StVKReducedInternalForces.h"
#include "volumetricMeshENuMaterial.h"

StVKReducedInternalForces::StVKReducedInternalForces(int r, double * U, VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, int initOnly, bool addGravity_, double g_, int verbose_): precomputedIntegrals(precomputedABCDIntegrals), unitReducedGravityForce(NULL), reducedGravityForce(NULL), addGravity(addGravity_), g(g_), useSingleThread(0), shallowCopy(0), verbose(verbose_)
{
  int numElements = volumetricMesh->getNumElements();
  lambdaLame = (double*) malloc (sizeof(double) * numElements);
  muLame = (double*) malloc (sizeof(double) * numElements);

  for(int el=0; el<numElements; el++)
  {
    VolumetricMesh::Material * material = volumetricMesh->getElementMaterial(el);
    VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
    if (eNuMaterial == NULL)
    {
      printf("Error: StVKReducedInternalForces: mesh does not consist of E, nu materials.\n");
      throw 1;
    }

    lambdaLame[el] = eNuMaterial->getLambda();
    muLame[el] = eNuMaterial->getMu();
  }

  InitComputation(r, U, volumetricMesh);
  if (!initOnly)
    ProcessElements(0, volumetricMesh->getNumElements());
  InitGravity();
}

StVKReducedInternalForces::StVKReducedInternalForces(const char * filename, int rTarget, int bigEndianMachine, int verbose_) : verbose(verbose_)
{
  FILE * fin = fopen(filename, "rb");
  if (!fin)
  {
    printf("Error: could not read from the input cubic polynomial file.\n");
    throw 1;
  }

  LoadFromStream(fin, rTarget, bigEndianMachine);
  fclose(fin);
}

StVKReducedInternalForces::StVKReducedInternalForces(FILE * fin, int rTarget, int bigEndianMachine, int verbose_) : verbose(verbose_)
{
  LoadFromStream(fin, rTarget, bigEndianMachine); 
}

int StVKReducedInternalForces::LoadFromStream(FILE * fin, int rTarget, int bigEndianMachine) 
{
  if (verbose)
    printf("Loading polynomials assuming little endian machine: %s.", (!bigEndianMachine) ? "TRUE" : "FALSE");

  int header[4];

  if ((int)(fread(header, sizeof(int), 4, fin)) < 4)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    throw 1;
  }
  
  r = header[0];

  int buffer;
  if (bigEndianMachine)
  {
    little2big(&r, &buffer, sizeof(int));
    r = buffer;
  }

  if (rTarget > r)
  {
    printf("Error: the input cubic polynomial file has r=%d, but you requested %d > %d.\n", r, rTarget, r);
    throw 2;
  }

  // first read in the coefficients as if all modes requested
  if (verbose)
    printf(" r=%d\n", r);
  
  r2 = r * r;

  linearSize = header[1];

  if (bigEndianMachine)
  {
    little2big(&linearSize, &buffer, sizeof(int));
    linearSize = buffer;
  }

  quadraticSize = header[2];

  if (bigEndianMachine)
  {
    little2big(&quadraticSize, &buffer, sizeof(int));
    quadraticSize = buffer;
  }

  cubicSize = header[3];

  if (bigEndianMachine)
  {
    little2big(&cubicSize, &buffer, sizeof(int));
    cubicSize = buffer;
  }

  linearCoef_ = (double*) malloc (sizeof(double) * r * linearSize);

  if ((int)(fread(linearCoef_,sizeof(double),r*linearSize,fin)) < r*linearSize)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    throw 1;
  }

  double bufferd;
  if (bigEndianMachine)
  {
    for(int i=0; i<r*linearSize; i++)
    {
      little2big(&linearCoef_[i], &bufferd, sizeof(double));
      linearCoef_[i] = bufferd;
    }
  }

  quadraticCoef_ = (double*) malloc (sizeof(double) * r * quadraticSize);

  if ((int)(fread(quadraticCoef_,sizeof(double),r*quadraticSize,fin)) < r*quadraticSize)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    throw 1;
  }

  if (bigEndianMachine)
  {
    for(int i=0; i<r*quadraticSize; i++)
    {
      little2big(&quadraticCoef_[i], &bufferd, sizeof(double));
      quadraticCoef_[i] = bufferd;
    }
  }

  cubicCoef_ = (double*) malloc (sizeof(double) * r * cubicSize);

  if ((int)(fread(cubicCoef_,sizeof(double),r*cubicSize,fin)) < r*cubicSize)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    throw 1;
  }

  if (bigEndianMachine)
  {
    for(int i=0; i<r*cubicSize; i++)
    {
      little2big(&cubicCoef_[i], &bufferd, sizeof(double));
      cubicCoef_[i] = bufferd;
    }
  }

  if (rTarget >= 0)
  {
    int linearSizeTarget, quadraticSizeTarget, cubicSizeTarget;
    GetSizes(rTarget, &linearSizeTarget, &quadraticSizeTarget, &cubicSizeTarget);

    double * linearCoefTemp_ = 
      (double*) malloc (sizeof(double) * rTarget * linearSizeTarget);

    double * quadraticCoefTemp_ = 
      (double*) malloc (sizeof(double) * rTarget * quadraticSizeTarget);

    double * cubicCoefTemp_ = 
      (double*) malloc (sizeof(double) * rTarget * cubicSizeTarget);

    for(int output=0; output<rTarget; output++)
      for(int i=0; i<rTarget; i++)
      {
        SetSizes(rTarget);
        int positionTarget = linearCoefPos(output, i); 
        SetSizes(r);
        int position = linearCoefPos(output, i); 
        linearCoefTemp_[positionTarget] = linearCoef_[position];
      }
 
    for(int output=0; output<rTarget; output++)
      for(int i=0; i<rTarget; i++)
        for(int j=i; j<rTarget; j++)
        {
          SetSizes(rTarget);
          int positionTarget = quadraticCoefPos(output, i, j); 
          SetSizes(r);
          int position = quadraticCoefPos(output, i, j); 
          quadraticCoefTemp_[positionTarget] = quadraticCoef_[position];
        }

    for(int output=0; output<rTarget; output++)
      for(int i=0; i<rTarget; i++)
        for(int j=i; j<rTarget; j++)
          for(int k=j; k<rTarget; k++)
          {
            SetSizes(rTarget);
            int positionTarget = cubicCoefPos(output, i, j, k); 
            SetSizes(r);
            int position = cubicCoefPos(output, i, j, k); 
            cubicCoefTemp_[positionTarget] = cubicCoef_[position];
          }

    r = rTarget;
    SetSizes(r);

    free(linearCoef_);
    free(quadraticCoef_);
    free(cubicCoef_);

    linearCoef_ = linearCoefTemp_;
    quadraticCoef_ = quadraticCoefTemp_;
    cubicCoef_ = cubicCoefTemp_;
  }

  volumetricMesh = NULL;
  U = NULL;
  reducedGravityForce = NULL;
  precomputedIntegrals = NULL;
  numElementVertices = 0;

  InitBuffers();

  addGravity = false;

  useSingleThread = 0;
  shallowCopy = 0;
  g=9.81; 

  return 0;
}

StVKReducedInternalForces::~StVKReducedInternalForces()
{
  if (!shallowCopy)
  {
    free(unitReducedGravityForce);
    free(reducedGravityForce);
    free(linearCoef_);
    free(quadraticCoef_);
    free(cubicCoef_);
    free(lambdaLame);
    free(muLame);
  }
  FreeBuffers();
}

void StVKReducedInternalForces::InitGravity(VolumetricMesh * volumetricMesh_, double * U_)
{
  VolumetricMesh * mesh = volumetricMesh;
  if (volumetricMesh_ != NULL)
    mesh = volumetricMesh_;
  
  double * UB = NULL;
  if (U != NULL)
    UB = U;
  if (U_ != NULL)
    UB = U_;

  if ((mesh == NULL) || (UB ==NULL))
  {
    printf("Error: cannot init gravity. Mesh or basis is not specified.\n");
    exit(1);
  }

  if (reducedGravityForce == NULL)
  {
    int n = mesh->getNumVertices();
    double * gravityForce = (double*) malloc (sizeof(double) * 3 * n);
    unitReducedGravityForce = (double*) malloc (sizeof(double) * r);
    reducedGravityForce = (double*) malloc (sizeof(double) * r);
    mesh->computeGravity(gravityForce, 1.0);
    ProjectVector(3*n, r, UB, unitReducedGravityForce, gravityForce);
    //for(int i=0; i<r; i++)
      //printf("%G\n", unitReducedGravityForce[i]);
    free(gravityForce);
  }

  for(int i=0; i<r; i++)
    reducedGravityForce[i] = g * unitReducedGravityForce[i];

  //printf("Altered gravity to %G\n", g);
}

void StVKReducedInternalForces::InitComputation(int r, double * U, VolumetricMesh * volumetricMesh)
{
  this->volumetricMesh = volumetricMesh;
  this->U = U;
  this->r = r;
  numElementVertices = volumetricMesh->getNumElementVertices();

  n = volumetricMesh->getNumVertices();

  r2 = r * r;

  // allocate room for coefficients, r linear coefficients per each of the r components
  linearSize = GetLinearSize(r);
  linearCoef_ = (double*) calloc (r * linearSize, sizeof(double));

  // allocate room for coefficients, r*(r+1)/2 quadratic coefficients per each of the r components
  quadraticSize = GetQuadraticSize(r);
  quadraticCoef_ = (double*) calloc ( r * quadraticSize, sizeof(double));

  // allocate room for coefficients, r*(r+1)*(r+2)/6 cubic coefficients per each of the r components
  cubicSize = GetCubicSize(r);
  cubicCoef_ = (double*) calloc (r * cubicSize, sizeof(double));

  InitBuffers();

  if (verbose >= 1)
  {
    printf("Number of vertices is: %d\n",n);
    printf("Number of nonlinear modes is: %d\n",r);
    printf("Computation initialization completed.\n");
  }
}

void StVKReducedInternalForces::GetSizes(int r, int * linearSize, int * quadraticSize, int * cubicSize)
{
  *linearSize = r;
  *quadraticSize = r * (r+1) / 2;
  *cubicSize = r * (r+1) * (r+2) / 6;
}

void StVKReducedInternalForces::ProcessElements(int startElement, int endElement, double ** target)
{
  double * linearCoef_ = this->linearCoef_;
  double * quadraticCoef_ = this->quadraticCoef_;
  double * cubicCoef_ = this->cubicCoef_;

  if (target != NULL)
  {
    linearCoef_ = target[0];
    quadraticCoef_ = target[1];
    cubicCoef_ = target[2];
  }

  if (verbose >= 1)
    printf("Generating element data: element %d to %d...\n", startElement, endElement-1);

  int numVertices_ = volumetricMesh->getNumVertices();

  // make auxiliary vectors
  double * qiqjBuffer = (double*) calloc(r2,sizeof(double));
  double * qkBuffer = (double*) calloc(r2,sizeof(double));
  double * coefs = (double*) calloc(r*r*r*r,sizeof(double));

  void * elIter;
  precomputedIntegrals->AllocateElementIterator(&elIter);

  // Linear terms
  //if (verbose >= 1)
    //printf("Building linear terms:");

  for(int el=startElement; el < endElement; el++)
  {
    precomputedIntegrals->PrepareElement(el, elIter);

    if (verbose >= 1)
    {
      if (el % 100 == 1)
        printf("%d ",el); fflush(NULL);
    }

    double lambda = lambdaLame[el];
    double mu = muLame[el];

    for(int i=0; i<r; i++)
    {
      for (int c=0; c<numElementVertices; c++)
      {
        Vec3d force(0.0,0.0,0.0);

        int vc = volumetricMesh->getVertexIndex(el, c);
        for (int a=0; a<numElementVertices; a++)
        {
          int va = volumetricMesh->getVertexIndex(el, a);

          Vec3d ua(U[ELT(3*numVertices_,3*va+0,i)],
                   U[ELT(3*numVertices_,3*va+1,i)],
                   U[ELT(3*numVertices_,3*va+2,i)]);

          force += lambda * (precomputedIntegrals->A(elIter,c,a) * ua) +
                   (mu * precomputedIntegrals->B(elIter,a,c)) * ua +
                   mu * (precomputedIntegrals->A(elIter,a,c) * ua);
        }

        // multiply Uc^T * force
        for(int output=0; output<r; output++)
        {
          linearCoef_[linearCoefPos(output, i)] +=
                  U[ELT(3*numVertices_,3*vc+0,output)] * force[0] +
                  U[ELT(3*numVertices_,3*vc+1,output)] * force[1] +
                  U[ELT(3*numVertices_,3*vc+2,output)] * force[2];
        }
      }
    }
  }

  // Quadratic terms
  //if (verbose >= 1)
    //printf("\nBuilding quadratic terms:");

  double ** forceBuffer = (double**) malloc (sizeof(double*) * numElementVertices);
  for(int c=0; c<numElementVertices; c++)
    forceBuffer[c] = (double*) calloc (3*r2,sizeof(double));

  memset(quadraticCoef_, 0, sizeof(double) * r * quadraticSize);

  int * vertices = (int*) malloc (sizeof(int) * numElementVertices);

  for(int el=startElement; el < endElement; el++)
  {
    precomputedIntegrals->PrepareElement(el, elIter);

    if (verbose >= 1)
    {
      if (el % 100 == 1)
        printf("%d ",el); fflush(NULL);
    }

    double lambda = lambdaLame[el];
    double mu = muLame[el];

    for(int ver=0; ver<numElementVertices ;ver++)
      vertices[ver] = volumetricMesh->getVertexIndex(el, ver);

    for(int c=0; c<numElementVertices; c++)
      memset(forceBuffer[c],0,sizeof(double)*3*r2);

    for(int a=0; a<numElementVertices; a++)
    {
      for(int b=0; b<numElementVertices; b++)
      {
        // compute ua*ub for all possible i,j
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                      r, r, 3,
                      1.0,
                      &U[ELT(3*numVertices_,3*vertices[a],0)], 3*numVertices_,
                      &U[ELT(3*numVertices_,3*vertices[b],0)], 3*numVertices_,
                      0.0,
                      qiqjBuffer, r);

        for(int c=0; c<numElementVertices; c++)
        {
          Vec3d vec1 = 0.5 * lambda * precomputedIntegrals->C(elIter,c,a,b) +
                       mu * precomputedIntegrals->C(elIter,a,b,c);

          Vec3d C = lambda * precomputedIntegrals->C(elIter,a,b,c) +
                    mu * (precomputedIntegrals->C(elIter,c,a,b) + precomputedIntegrals->C(elIter,b,a,c)); 

          for(int i=0; i<r; i++)
          {
            double * posa = &(U[ELT(3*numVertices_,3*vertices[a]+0,i)]);
            double Cdotua = C[0] * posa[0] + C[1] * posa[1] + C[2] * posa[2];

            for(int j=0; j<r; j++)
            {
              double buffer = qiqjBuffer[ELT(r,i,j)];
              double * posb = &(U[ELT(3*numVertices_,3*vertices[b]+0,j)]);

              int index = ELT(3,0,ELT(r,i,j));

              forceBuffer[c][index+0] += buffer * vec1[0] + Cdotua * posb[0];
              forceBuffer[c][index+1] += buffer * vec1[1] + Cdotua * posb[1];
              forceBuffer[c][index+2] += buffer * vec1[2] + Cdotua * posb[2];
            }
          }
        } // end c
      } // end b
    } // end a

    // generate unpacked coefficients for this element
    memset(coefs,0,sizeof(double)*r*r*r);
    for(int c=0; c<numElementVertices; c++)
    {
      // multiply Uc^T * forcesBuffer[c]
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    r, r2, 3,
                    1.0,
                    &U[ELT(3*numVertices_,3*vertices[c],0)], 3*numVertices_,
                    forceBuffer[c], 3,
                    1.0,
                    coefs, r);
    }

    // pack and add
    for(int output=0; output<r; output++)
    {
      for(int i=0; i<r; i++)
        for(int j=0; j<r; j++)
        {
          int i1 = i;
          int j1 = j;

          if (j < i)
          {
            i1 = j;
            j1 = i;
          }

          quadraticCoef_[quadraticCoefPos(output,i1, j1)] += coefs[ELT(r,output,ELT(r,i,j))];
        }
    }
  } // end el

  free(vertices);

  for(int c=0; c<numElementVertices; c++)
    free(forceBuffer[c]);
  free(forceBuffer);

  // cubic terms
  //if (verbose >= 1)
    //printf("\nBuilding cubic terms:\n");

  memset(coefs,0,sizeof(double)*r*r*r*r);

  for(int el=startElement; el < endElement; el++)
  {
    precomputedIntegrals->PrepareElement(el, elIter);

    if (verbose >= 1)
    {
      if ((el % 50 == 1) || ((r > 30) && (el % 25 == 1)))
        printf("%d ",el); fflush(NULL);
    }

    double lambda = lambdaLame[el];
    double mu = muLame[el];

    for (int a=0; a<numElementVertices; a++)
    {
      int va = volumetricMesh->getVertexIndex(el, a);
      for(int b=0; b<numElementVertices; b++)
      {
        int vb = volumetricMesh->getVertexIndex(el, b);

        // fill up the buffers
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                         r, r, 3,
                         1.0,
                         &U[ELT(3*numVertices_,3*va,0)], 3*numVertices_,
                         &U[ELT(3*numVertices_,3*vb,0)], 3*numVertices_,
                         0.0,
                         qiqjBuffer, r);

        for(int i=0; i<r2; i++)
          qkBuffer[i] = 0;

        for(int c=0; c<numElementVertices; c++)
        {
          int vc = volumetricMesh->getVertexIndex(el, c);
          for(int d=0; d<numElementVertices; d++)
          {
            int vd = volumetricMesh->getVertexIndex(el, d);

            double factor = 0.5 * lambda * precomputedIntegrals->D(elIter,a,b,c,d) +
                            mu * precomputedIntegrals->D(elIter,a,c,b,d);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            r, r, 3,
                            factor,
                            &U[ELT(3*numVertices_,3*vc,0)], 3*numVertices_,
                            &U[ELT(3*numVertices_,3*vd,0)], 3*numVertices_,
                            1.0,
                            qkBuffer, r);
          }
        }

        // multiply qiqjBuffer * qkBuffer^T (tensor product)
        // both vectors are r^2 vectors

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                          r2, r2, 1,
                          1.0,
                          qiqjBuffer, r2,
                          qkBuffer, r2,
                          1.0,
                          coefs, r2);

      } // over b
    } // over a
  }

  for(int i=0; i < r*cubicSize; i++)
    cubicCoef_[i] = 0.0;

  // unpack
  for(int i=0; i<r; i++)
    for(int j=0; j<r; j++)
      for(int k=0; k<r; k++)
        for(int l=0; l<r; l++)
        {
          // sort the indices
          int i1=i;
          int j1=j;
          int k1=k;
          tripleSort(i1,j1,k1);

          //int pos = cubicCoefPos(l, i1, j1, k1);
          //int pos1 = ELT(r*r,ELT(r,i,j),ELT(r,l,k));
          cubicCoef_[cubicCoefPos(l, i1, j1, k1)] += coefs[ELT(r*r,ELT(r,i,j),ELT(r,l,k))];
        }

  free(qiqjBuffer);
  free(qkBuffer);
  free(coefs);

  precomputedIntegrals->ReleaseElementIterator(elIter);
}

void StVKReducedInternalForces::TestPolynomial(double * q, StVKInternalForces * stVKInternalForces, const char * filenameU)
{
  double * fqPoly = (double*) malloc (sizeof(double) * r);
  double * fqDirect = (double*) malloc (sizeof(double) * r);

  Evaluate(q, fqPoly);

  double * U1;
  int n1,r1;

  if (ReadMatrixFromDisk(filenameU, &n1, &r1, &U1) != 0)
  {
    printf("Error loading file %s\n",filenameU);
    exit(1);
  }

  n1 /= 3;

  printf("\nNumber of vertices is: %d\n",n1);
  printf("Number of nonlinear modes is: %d\n",r1);

  double * u = (double*) malloc (sizeof(double) * 3 * n1);
  double * forces = (double*) malloc (sizeof(double) * 3 * n1);

  // u = U*q
  SynthesizeVector(3*n1, r1, U1, q, u);

  stVKInternalForces->ComputeForces(u,forces);

/*
  printf("g***\n");
  for(int j=0; j<3*n1; j++)
    printf("%G\n",forces[j]);
  printf("g***\n");
*/

  ProjectVector(3*n1,r1,U1,fqDirect,forces);

  int i;

  printf("ViaCubPoly ViaProjection:\n");

  for(i=0; i<r1; i++)
  {
    printf("%G %G\n",fqPoly[i],fqDirect[i]);
  }

  free(u);
  free(forces);

  free(U1);

  free(fqDirect);
  free(fqPoly);
}

int StVKReducedInternalForces::Save(const char * filename)
{
  FILE * fout = fopen(filename,"wb");
  if (!fout)
    return 1;

  Save(fout);

  fclose(fout);

  return 0;
}

int StVKReducedInternalForces::Save(FILE * fout)
{
  if ((int)(fwrite(&r,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(&linearSize,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(&quadraticSize,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(&cubicSize,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(linearCoef_,sizeof(double),r*linearSize,fout)) < r*linearSize)
    return 1;

  if ((int)(fwrite(quadraticCoef_,sizeof(double),r*quadraticSize,fout)) < r*quadraticSize)
    return 1;

  if ((int)(fwrite(cubicCoef_,sizeof(double),r*cubicSize,fout)) < r*cubicSize)
    return 1;

  return 0;
}

void StVKReducedInternalForces::SetSizes(int rTarget)
{
  r = rTarget;
  r2 = r*r;
  GetSizes(r, &linearSize, &quadraticSize, &cubicSize);
}

int StVKReducedInternalForces::GetrFromFile(const char * filename)
{
  FILE * fin = fopen(filename,"rb");

  if (!fin)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    return -1; 
  }

  int r;
  if ((int)(fread(&r,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input cubic polynomial file.\n");
    return -1;
  }

  fclose(fin);

  return r;
}

void StVKReducedInternalForces::little2big(void * input, void * output, int numBytes)
{
  char * inp = (char*) input;
  char * out = (char*) output;
  out += numBytes;
  for(int i=0; i<numBytes; i++)
  {
    out--;
    *out = *inp;
    inp++;
  }
}


void StVKReducedInternalForces::Evaluate(double * q, double * fq)
{
/* // unoptimized version
  // reset to zero
  int i,j,k,l;
  for(l=0; l<r; l++)
    fq[l] = 0;

  // add linear terms
  int index = 0;
  for(l=0; l<r; l++)
    for(i=0; i<r; i++)
    {
      fq[l] += linearCoef_[index] * q[i];
      index++;
    }

  // add quadratic terms
  index = 0;
  for(l=0; l<r; l++)
    for(i=0; i<r; i++)
      for(j=i; j<r; j++)
      {
        fq[l] += quadraticCoef_[index] * q[i] * q[j];
        index++;
      }

  // add cubic terms
  index = 0;
  for(l=0; l<r; l++)
    for(i=0; i<r; i++)
      for(j=i; j<r; j++)
        for(k=j; k<r; k++)
        {
          fq[l] += cubicCoef_[index] * q[i] * q[j] * q[k];
          index++;
        }
*/

  if (useSingleThread)
  {
    #if defined(WIN32) || defined(linux)
      mkl_max_threads = mkl_get_max_threads();
      mkl_dynamic = mkl_get_dynamic();
      mkl_set_num_threads(1);
      mkl_set_dynamic(0);
    #elif defined(__APPLE__)
      //setenv("VECLIB_MAXIMUM_THREADS", "1", true);
    #endif
  }

  // add linear terms
  // multiply linearCoef_ and q
  // linearCoef_ is r x r array
  cblas_dgemv(CblasColMajor, CblasTrans,
       r, r,
       1.0,
       linearCoef_, r,
       q, 1,
       0.0,
       fq, 1);

  // compute qiqj
  int index = 0;
  for(int output=0; output<r; output++)
    for(int i=output; i<r; i++)
    {
      qiqj[index] = q[output] * q[i];
      index++;
    }

  // add quadratic terms
  // quadraticCoef_ is quadraticSize x r matrix
  // each column gives quadratic coef for one force vector component
  cblas_dgemv(CblasColMajor, CblasTrans,
       quadraticSize, r,
       1.0,
       quadraticCoef_, quadraticSize,
       qiqj, 1,
       1.0,
       fq, 1);

  // add cubic terms
  // cubicCoef_ is cubicSize x r matrix
  // each column gives cubicSize coef for one force vector component
  int size = quadraticSize;
  double * qiqjPos = qiqj;
  double * cubicCoefPos = cubicCoef_;
  for(int i=0; i<r; i++)
  {
    cblas_dgemv(CblasColMajor, CblasTrans,
        size, r,
        q[i],
        cubicCoefPos, cubicSize,
        qiqjPos, 1,
        1.0,
        fq, 1);

    int param = r-i;
    size -= param;
    qiqjPos += param;
    cubicCoefPos += param * (param+1) / 2;
  }

  if (addGravity)
  {
    for(int i=0; i<r; i++)
      fq[i] -= reducedGravityForce[i];
  }

  if (useSingleThread)
  {
    #if defined(WIN32) || defined(linux)
      mkl_set_num_threads(mkl_max_threads);
      mkl_set_dynamic(mkl_dynamic);
    #elif defined(__APPLE__)
      //unsetenv("VECLIB_MAXIMUM_THREADS");
    #endif
  }
}

double StVKReducedInternalForces::EvaluateComponent(double * q, int componentIndex)
{
  // add linear terms
  // multiply linearCoef_ and q
  // linearCoef_ is r x r array

  double result = cblas_ddot(r, linearCoef_ + componentIndex * r, 1, q, 1);

  // compute qiqj
  int index = 0;
  for(int output=0; output<r; output++)
    for(int i=output; i<r; i++)
    {
      qiqj[index] = q[output] * q[i];
      index++;
    }

  // add quadratic terms
  // quadraticCoef_ is quadraticSize x r matrix
  // each column gives quadratic coef for one force vector component

  result += cblas_ddot(quadraticSize, 
    quadraticCoef_ + componentIndex * quadraticSize, 1,
    qiqj, 1);

  // add cubic terms
  // cubicCoef_ is cubicSize x r matrix
  // each column gives cubicSize coef for one force vector component
  int size = quadraticSize;
  double * qiqjPos = qiqj;
  double * cubicCoefPos = cubicCoef_ + cubicSize * componentIndex;
  for(int i=0; i<r; i++)
  {
    result += q[i] * cblas_ddot(size, cubicCoefPos, 1, qiqjPos, 1);
    int param = r-i;
    size -= param;
    qiqjPos += param;
    cubicCoefPos += param * (param+1) / 2;
  }

  return result;
}

void StVKReducedInternalForces::EvaluateLinear(double * q, double * fq)
{
  // reset to zero
  for(int l=0; l<r; l++)
    fq[l] = 0;

  // add linear terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
    {
      fq[l] += linearCoef_[index] * q[i];
      index++;
    }
}

void StVKReducedInternalForces::PrintPolynomial()
{
  // linear terms
  int index = 0;
  printf("{");
  for(int l=0; l<r; l++)
  {
    for(int i=0; i<r; i++)
    {
      printf("%.15f * q%d",linearCoef_[index],i);
      if (i != r - 1)
        printf(" + ");
      index++;
      if (i % 5 == 4)
        printf("\n");
    }
    if (l != r-1)
      printf(", ");
  }
  printf("} + ");

  // quadratic terms
  index = 0;
  printf("{");
  for(int l=0; l<r; l++)
  {
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
      {
        printf("%.15f * q%d * q%d",quadraticCoef_[index],i,j);
        index++;
        if (!((i == r - 1) && (j == r - 1)))
          printf(" + ");

        if (index % 5 == 4)
          printf("\n");

      }
    if (l != r-1)
      printf(", ");
  }
  printf("} + ");

  // cubic terms
  index = 0;
  printf("{");
  for(int l=0; l<r; l++)
  {
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          printf("%.15f * q%d * q%d * q%d",cubicCoef_[index],i,j,k);

          if (!((i == r - 1) && (j == r - 1) && (k == r - 1)))
            printf(" + ");

          if (index % 5 == 4)
            printf("\n");

          index++;
        }

    if (l != r-1)
      printf(", ");
  }
  printf("}");
}

void StVKReducedInternalForces::tripleSort(int & i, int & j, int & k)
{
  // sort in ascending order

  int buffer;
  #define SWAP(a,b)\
  buffer = a;\
  a = b;\
  b = buffer;

  // bubble sort on 3 elements
  if (j < i)
  {
    SWAP(i,j);
  }

  if (k < j)
  {
    SWAP(j,k);
  }

  if (j < i)
  {
    SWAP(i,j);
  }


}

int StVKReducedInternalForces::CubicCoefficientsHistogramLog10(int * histogram,int lowExp, int highExp)
{
  int flag = 0;

  // set to zero
  for(int i=0; i<highExp-lowExp+1; i++)
    histogram[i] = 0;

  // over all cubic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          double term = fabs(cubicCoef_[index]);

          // determine bucket
          double log10term =  log10(term);
          int log10discrete = (int)log10term;

          if (log10discrete > highExp)
          {
            printf("Warning: encountered entry above given histogram interval.\n");
            flag = 1;
            continue;
          }

          if (log10discrete < lowExp)
          {
            printf("Warning: encountered entry below given histogram interval.\n");
            flag = 1;
            continue;
          }

          histogram[log10discrete-lowExp]++;

          index++;
        }

  return flag;

}

// finds max coefficient in terms of abs value among all terms
double StVKReducedInternalForces::MaxAbsCoefficient()
{
  double max = 0;

  // linear terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
    {
      double term = fabs(linearCoef_[index]);
      if (term > max)
        max = term;
      index++;
    }

  // quadratic terms
  index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
      {
        double term = fabs(quadraticCoef_[index]);
        if (term > max)
          max = term;
        index++;
      }

  // cubic terms
  index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          double term = fabs(cubicCoef_[index]);
          if (term > max)
            max = term;
          index++;
        }

  return max;
}

// finds max coefficient in terms of abs value among all terms
double StVKReducedInternalForces::MaxAbsLinearCoefficient()
{
  double max = 0;

  // over linear terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
    {
      double term = fabs(linearCoef_[index]);
      if (term > max)
        max = term;
      index++;
    }

  return max;
}

// finds max coefficient in terms of abs value among all terms
double StVKReducedInternalForces::MaxAbsQuadraticCoefficient()
{
  double max = 0;

  // over quadratic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
      {
        double term = fabs(quadraticCoef_[index]);
        if (term > max)
          max = term;
        index++;
      }

  return max;
}

// finds max coefficient in terms of abs value among all terms
double StVKReducedInternalForces::MaxAbsCubicCoefficient()
{
  double max = 0;

  // over cubic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          double term = fabs(cubicCoef_[index]);
          if (term > max)
            max = term;
          index++;
        }

  return max;
}

double StVKReducedInternalForces::AverageAbsCubicCoefficient()
{
  double avg = 0;

  // add cubic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          avg += fabs(cubicCoef_[index]);
          index++;
        }

  return (1.0 * avg / index);
}

void StVKReducedInternalForces::PrintLinearCoefficients()
{
  // linear terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
    {
      double term = linearCoef_[index];
      //term = fabs(term);
      printf("%G\n",term);
      index++;
    }
}

void StVKReducedInternalForces::PrintQuadraticCoefficients()
{
  // quadratic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
      {
        double term = quadraticCoef_[index];
        //term = fabs(term);
        printf("%G\n",term);
        index++;
      }
}

// finds max coefficient in terms of abs value among all terms
void StVKReducedInternalForces::PrintCubicCoefficients()
{
  // add cubic terms
  int index = 0;
  for(int l=0; l<r; l++)
    for(int i=0; i<r; i++)
      for(int j=i; j<r; j++)
        for(int k=j; k<r; k++)
        {
          double term = cubicCoef_[index];
          //term = fabs(term);
          printf("%G\n",term);
          index++;
        }
}

void StVKReducedInternalForces::UseSingleThread(int useSingleThread_)
{
  useSingleThread = useSingleThread_;
}

void StVKReducedInternalForces::InitBuffers()
{
  qiqj = (double*) malloc (sizeof(double) * quadraticSize);
}

void StVKReducedInternalForces::FreeBuffers()
{
  free(qiqj);
}

void StVKReducedInternalForces::Scale(double scalingFactor)
{
  for(int i=0; i<r*linearSize; i++)
    linearCoef_[i] *= scalingFactor;

  for(int i=0; i<r*quadraticSize; i++)
    quadraticCoef_[i] *= scalingFactor;

  for(int i=0; i<r*cubicSize; i++)
    cubicCoef_[i] *= scalingFactor;
}

StVKReducedInternalForces * StVKReducedInternalForces::ShallowClone()
{
  StVKReducedInternalForces * output = new StVKReducedInternalForces(*this); // invoke default copy constructor
  output->shallowCopy = 1;
  output->InitBuffers();
  return output;
}

int StVKReducedInternalForces::SaveEmptyCub(const char * filename)
{
  FILE * fout = fopen(filename, "wb");
  if (!fout)
    return 1;

  SaveEmptyCub(fout);

  fclose(fout);
  return 0;
}

int StVKReducedInternalForces::SaveEmptyCub(FILE * fout)
{
  int r = 0;
  if ((int)(fwrite(&r,sizeof(int),1,fout)) < 1)
    return 1;

  int linearSize = 0;
  if ((int)(fwrite(&linearSize,sizeof(int),1,fout)) < 1)
    return 1;

  int quadraticSize = 0;
  if ((int)(fwrite(&quadraticSize,sizeof(int),1,fout)) < 1)
    return 1;

  int cubicSize = 0;
  if ((int)(fwrite(&cubicSize,sizeof(int),1,fout)) < 1)
    return 1;

  return 0;
}

/*
  // the cached version of Evaluate... performs slowly than the straight version
  // build products q_i * q_j
  index = 0;
  for(i=0; i<r; i++)
    for(j=i; j<r; j++)
    {
      qij[index] = q[i] * q[j];
      index++;
    }


  // add quadratic terms
  index = 0;
  double * qijData;
  for(l=0; l<r; l++)
  {
    qijData = qij;
    for(i=0; i<r; i++)
      for(j=i; j<r; j++)
      {
        fq[l] += quadraticCoef_[index] * (*qijData);
        index++;
        qijData++;
      }
  }

  // add cubic terms
  index = 0;
  for(l=0; l<r; l++)
  {
    qijData = qij;
    for(i=0; i<r; i++)
      for(j=i; j<r; j++)
      {
        for(k=j; k<r; k++)
        {
          fq[l] += cubicCoef_[index] * (*qijData) * q[k];
          index++;
        }
        qijData++;
      }
  }
  */

