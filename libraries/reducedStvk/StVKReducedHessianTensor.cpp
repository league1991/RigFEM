/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
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
#include "matrixMacros.h"
#include "StVKReducedHessianTensor.h"

StVKReducedHessianTensor::StVKReducedHessianTensor(StVKReducedStiffnessMatrix * stVKReducedStiffnessMatrix) : shallowCopy(0)
{
  r = stVKReducedStiffnessMatrix->Getr();
  r2 = r*r;

  InitBuffers();

  printf("Building the reduced Hessian tensor linear polynomials... r is %d\n",r);

  int i;

  int output,input,deriv;

  linearSize = StVKReducedInternalForces::GetLinearSize(r);
  quadraticSize = StVKReducedInternalForces::GetQuadraticSize(r);

  printf("Building free terms:");

  // free terms
  // allocate room for coefficients
  freeCoef_ = (double*) malloc (sizeof(double) * quadraticSize * r);

  // obtain free terms by analytic derivation of the linear stiffness matrix terms
  for(output=0; output<r; output++)
  {
    printf(" %d",output);
    for(input=output; input<r; input++)
    {
      for(deriv=0; deriv<r; deriv++)
      {
        int pos = freeCoefPos(output,input,deriv);
        freeCoef_[pos] = stVKReducedStiffnessMatrix->linearCoef(output,input,deriv);
      }
    }
  }

  printf("\nBuilding linear terms:");

  // linear terms
  // allocate room for coefficients, r coefficients per each of cubicSize components
  linearCoef_ = (double*) malloc (sizeof(double) * quadraticSize * r * r);

  // obtain linear coefficients by analytic derivation of the quadratic stiffness matrix terms
  for(output=0; output<r; output++)
  {
    printf(" %d",output);
    for(input=output; input<r; input++)
      for(deriv=0; deriv<r; deriv++)
      {
        for(i=0; i<r; i++)
        {
          // (i1,j1) will be (i,deriv) sorted in ascending order
          int i1 = i;
          int j1 = deriv;
          if (j1 < i1) // swap them
          {
            j1 = i;
            i1 = deriv;
          }

          double value = stVKReducedStiffnessMatrix->quadraticCoef(output,input,i1,j1);

          if (i == deriv)
            value *= 2;

          int pos = linearCoefPos(output,input,deriv,i);

          linearCoef_[pos] = value;
        }
      }
  }

  printf("\n");
}

StVKReducedHessianTensor::~StVKReducedHessianTensor()
{
  if (!shallowCopy)
  {
    free(freeCoef_);
    free(linearCoef_);
  }
  FreeBuffers();
}

void StVKReducedHessianTensor::Evaluate(double * q, double * Hq)
{
  // this is same as EvaluateSubset with start=0, end=quadraticSize

/*
  int i;
  int output,input,deriv;
*/
  // reset to free terms
  
  memcpy(Hq,freeCoef_,quadraticSize*r*sizeof(double));

/*
  int index = 0;
  for(output=0; output<r; output++)
  {
    for(input=output; input<r; input++)
    {
      for(deriv=input; deriv<r; deriv++)
      {
        Hq[index] = freeCoef_[index];
        index++;
      }
     // consider incrementing a separate index? 
    }
  }
*/
  // add linear terms
/*
void cblas_dgemv(const enum CBLAS_ORDER order, const enum
CBLAS_TRANSPOSE TransA, const int M, const int N, const double
alpha, const double *A, const int lda, const double *X, const
int incX, const double beta, double *Y, const int incY);
*/

  // multiply linearCoef_ and q
  cblas_dgemv(CblasColMajor, CblasTrans, 
       r, quadraticSize*r,
       1.0,
       linearCoef_, r,
       q, 1,
       1.0,
       Hq, 1);

  /*

  index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(input=output; input<r; input++)
    {
      for(deriv=input; deriv<r; deriv++)
      {
        for(i=0; i<r; i++)
        {
          Hq[indexEntry] += linearCoef_[index] * q[i];
          index++;
        }
        indexEntry++;
      }
    }
    indexEntry += output + 1;
  }
*/
/*
  // make symmetric
  indexEntry = 0;
  for(output=0; output<r; output++)
    for(i=0; i<output; i++)
      Rq[output * r + i] = Rq[i * r + output];
*/
}

// evaluates only entries [start,...end), in the upper-triangular order 
//   (0 <= start <= end <= quadraticSize)
void StVKReducedHessianTensor::EvaluateSubset(double * q, int start, int end, double * Hq)
{
/*
  int i;
  int output,input,deriv;
*/
  // reset to free terms
  memcpy(Hq + r * start,freeCoef_ + r * start, (end-start)*r*sizeof(double));

/*
  int index = 0;
  for(output=0; output<r; output++)
  {
    for(input=output; input<r; input++)
    {
      for(deriv=input; deriv<r; deriv++)
      {
        Hq[index] = freeCoef_[index];
        index++;
      }
      // incrementing a separate index? 
    }
  }
*/
  // add linear terms
/*
void cblas_dgemv(const enum CBLAS_ORDER order, const enum
CBLAS_TRANSPOSE TransA, const int M, const int N, const double
alpha, const double *A, const int lda, const double *X, const
int incX, const double beta, double *Y, const int incY);
*/

  // multiply linearCoef_ and q
  cblas_dgemv(CblasColMajor, CblasTrans, 
       r, (end-start)*r,
       1.0,
       linearCoef_ + r*r* start, r,
       q, 1,
       1.0,
       Hq + r * start, 1);

  /*
  index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(input=output; input<r; input++)
    {
      for(deriv=input; deriv<r; deriv++)
      {
        for(i=0; i<r; i++)
        {
          Hq[indexEntry] += linearCoef_[index] * q[i];
          index++;
        }
        indexEntry++;
      }
    }
    //indexEntry += output + 1;
  }
*/
/*
  // make symmetric
  indexEntry = 0;
  for(output=0; output<r; output++)
    for(i=0; i<output; i++)
      Rq[output * r + i] = Rq[i * r + output];
*/
}

void StVKReducedHessianTensor::ContractWithVector(int r, double * Hq, double * q, double * A)
{
  // computes A = Hq : q

  int quadraticSize = StVKReducedInternalForces::GetQuadraticSize(r);

  // multiply Hq and q
  cblas_dgemv(CblasColMajor, CblasTrans, 
       r, quadraticSize, 
       1.0,
       Hq, r,
       q, 1,
       0.0,
       A, 1);

  for(int j=r-1; j>=0; j--)
    for(int i=r-1; i>=j; i--)
    {
      int lowerTrianglePos = j * r - (j-1) * j / 2 + (i-j);
      A[ELT(r,i,j)] = A[lowerTrianglePos];
      A[ELT(r,j,i)] = A[lowerTrianglePos];
    }
}

void StVKReducedHessianTensor::MakeRoomForTensor(double ** Hq)
{
  *Hq = (double*) malloc (sizeof(double) * quadraticSize * r);
}

void StVKReducedHessianTensor::ApproximateReducedForceAndStiffnessMatrix
  (double * baseHq, double * baseRq, double * dq, double * dfq, double * dRq)
{
  /*
  memset(dfq,0,sizeof(double)*r);
  memset(dRq,0,sizeof(double)*r*r);
  return;
  */

  ContractWithVector(r, baseHq, dq, dRq);

/*
  // multiply Hq and dq
  cblas_dgemv(CblasColMajor, CblasTrans, 
       r, quadraticSize, 
       1.0,
       baseHq, r,
       dq, 1,
       0.0,
       buffer1, 1);

  // triangular buffer1 now contains upper triangle of dRq = H:dq

  // write buffer1 into r x r form (into dRq)
  int i1=0,j1=0;
  for(int i=0; i< quadraticSize; i++)
  {
    dRq[ELT(r,i1,j1)] = buffer1[i];
    dRq[ELT(r,j1,i1)] = buffer1[i];
    j1++;
    if(j1 == r)
    {
      i1++;
      j1 = i1;
    }
  }
*/

  // dRq now contains correct values

  // copy baseRq + 1/2 * H:dq into rxr buffer2
  /*
  void cblas_daxpy(const int N, const double alpha, const double
*X, const int incX, double *Y, const int incY);
*/
/*
void cblas_dcopy(const int N, const double *X, const int incX,
double *Y, const int incY);  
*/
  cblas_dcopy(r2, baseRq, 1, buffer2, 1);
  cblas_daxpy(r2, 0.5, dRq, 1, buffer2, 1);

  // multiply buffer2 * dq = dfq
  cblas_dgemv(CblasColMajor, CblasNoTrans, 
      r, r, 
      1.0,
      buffer2, r,
      dq, 1,
      0.0,
      dfq, 1);
}

int StVKReducedHessianTensor::Save(const char * filename)
{
  FILE * fout = fopen(filename,"wb");

  if (!fout)
    return 1;

  if ((int)(fwrite(&r,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(&linearSize,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(&quadraticSize,sizeof(int),1,fout)) < 1)
    return 1;

  if ((int)(fwrite(freeCoef_,sizeof(double),quadraticSize*r,fout)) < quadraticSize*r)
    return 1;

  if ((int)(fwrite(linearCoef_,sizeof(double),quadraticSize*r*r,fout)) < quadraticSize*r*r)
    return 1;

  fclose(fout);

  return 0;
}


StVKReducedHessianTensor::StVKReducedHessianTensor(const char * filename)
{
  FILE * fin = fopen(filename,"rb");

  if (!fin)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  if ((int)(fread(&r,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  r2 = r * r;

  if ((int)(fread(&linearSize,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  if ((int)(fread(&quadraticSize,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  freeCoef_ = (double*) malloc (sizeof(double) * quadraticSize*r);

  if ((int)(fread(freeCoef_,sizeof(double),quadraticSize*r,fin)) < quadraticSize*r)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  linearCoef_ = (double*) malloc (sizeof(double) * quadraticSize * r * r);

  if ((int)(fread(linearCoef_,sizeof(double),quadraticSize*r*r,fin)) < quadraticSize*r*r)
  {
    printf("Error: couldn't read from input Hessian tensor file.\n");
    return;
  }

  fclose(fin);

  InitBuffers();
}

void StVKReducedHessianTensor::PrintTensor()
{
  int i;
  int output;
  int input,deriv;

  // free terms
  int index = 0;
  printf("{");
  for(output=0; output<r; output++)
  {
    printf("{");
    for(input=0; input<r; input++)
    {
      int i1 = output;
      int j1 = input;
      if (j1 < i1) // swap them
      {
        j1 = output;
        i1 = input;
      }
      // now i1 <= j1

      printf("{");
      for(deriv=0; deriv<r; deriv++)
      {
        printf("%.15f",freeCoef(i1,j1,deriv));
        if (deriv != r-1)
          printf(", ");
        index++;
        if (index % 5 == 4)
          printf("\n");
      }
      printf("}");
      if (input != r-1)
        printf(", ");
    }
    printf("}");
    if (output != r-1)
      printf(", ");
  }

  printf("} +\n{");

  index = 0;
  for(output=0; output<r; output++)
  {
    printf("{");
    for(input=0; input<r; input++)
    {
      printf("{");
      int i1 = output;
      int j1 = input;
      if (j1 < i1) // swap them
      {
        j1 = output;
        i1 = input;
      }
      // now i1 <= j1

      for(deriv=0; deriv<r; deriv++)
      {
        for(i=0; i<r; i++)
        {
          printf("%.15f * q%d",linearCoef(i1,j1,deriv,i),i);
          if (i != r-1)
            printf("+ ");
          index++;
          if (index % 5 == 4)
            printf("\n");
        }
        if (deriv != r-1)
          printf(", ");
      }
      printf("}");
      if (input != r-1)
        printf(", ");
    }
    printf("}");
    if (output != r-1)
      printf(", ");
  }
  printf("}");
}

void StVKReducedHessianTensor::InitBuffers()
{
  buffer1 = (double*) calloc (r2,sizeof(double));
  buffer2 = (double*) calloc (r2,sizeof(double));
}

void StVKReducedHessianTensor::FreeBuffers()
{
  free(buffer1);
  free(buffer2);
}

StVKReducedHessianTensor * StVKReducedHessianTensor::ShallowClone()
{
  StVKReducedHessianTensor * output = new StVKReducedHessianTensor(*this); // invoke default copy constructor
  output->shallowCopy = 1;
  output->InitBuffers();
  return output;
}

