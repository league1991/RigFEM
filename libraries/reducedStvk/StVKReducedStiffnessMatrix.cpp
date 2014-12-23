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
#if defined(WIN32) || defined(linux)
  #include "mkl_service.h"
#endif
#include "matrixMacros.h"
#include "StVKReducedStiffnessMatrix.h"

StVKReducedStiffnessMatrix::~StVKReducedStiffnessMatrix()
{
  if (!shallowCopy)
  {
    free(freeCoef_);
    free(linearCoef_);
    free(quadraticCoef_);
  }
  FreeBuffers();
}

StVKReducedStiffnessMatrix::StVKReducedStiffnessMatrix(StVKReducedInternalForces * stVKReducedInternalForces, int verbose) : shallowCopy(0)
{
  r = stVKReducedInternalForces->Getr();
  r2 = r*r;

  if (verbose)
    printf("Building the reduced stiffness matrix quadratic polynomials... r is %d\n",r);

  int i,j,k;

  int output;

  if (verbose)
    printf("Building free terms:");

  // free terms
  // allocate room for coefficients, 1 coefficient per each of the r x r components
  freeCoef_ = (double*) malloc (sizeof(double) * r * (r+1) / 2);

  // obtain free terms by analytic derivation of the linear force terms
  for(output=0; output<r; output++)
  {
    if (verbose)
      printf(" %d",output);
    for(i=output; i<r; i++)
    {
      freeCoef_[freeCoefPos(output,i)] = stVKReducedInternalForces->linearCoef(output,i);
    }
  }

  if (verbose)
    printf("\nBuilding linear terms:");

  // linear terms
  // allocate room for coefficients, r coefficients per each of the r x r components
  linearSize = StVKReducedInternalForces::GetLinearSize(r);
  linearCoef_ = (double*) malloc (sizeof(double) * r * (r+1) / 2 * linearSize);

  // obtain linear coefficients by analytic derivation of the quadratic force terms
  for(output=0; output<r; output++)
  {
    if (verbose)
      printf(" %d",output);
    for(i=output; i<r; i++)
      for(j=0; j<r; j++)
      {
        // (i1,j1) will be (i,j) sorted in ascending order
        int i1 = i;
        int j1 = j;
        if (j1 < i1) // swap them
        {
          j1 = i;
          i1 = j;
        }

        double value = stVKReducedInternalForces->quadraticCoef(output,i1,j1);

        if (i == j)
          value *= 2;

        //int pos = linearCoefPos(output,i,j);
        linearCoef_[linearCoefPos(output,i,j)] = value;
    }
  }

  if (verbose)
    printf("\nBuilding quadratic terms:");

  // quadratic terms
  // allocate room for coefficients, r*(r+1)/2 coefficients per each of the r x r components
  quadraticSize = StVKReducedInternalForces::GetQuadraticSize(r);
  quadraticCoef_ = (double*) malloc (sizeof(double) * r * (r+1) / 2 * quadraticSize);

  // obtain quadratic coefficients by analytic derivation of the cubic force terms
  for(output=0; output<r; output++)
  {
    if (verbose)
      printf(" %d",output);

    for(i=output; i<r; i++)
      for(j=0; j<r; j++)
        for(k=j; k<r; k++)
        {
          // (i1,j1,k1) will be (i,j,k) sorted in ascending order

          int i1 = i;
          int j1 = j;
          int k1 = k;

          int buffer;
          #define SWAP(i,j)\
             buffer = i;\
             i = j;\
             j = buffer;

          // bubble sort on 3 elements
          if (j1 < i1) 
          {
            SWAP(i1,j1);
          }

          if (k1 < j1) 
          {
            SWAP(j1,k1);
          }

          if (j1 < i1)
          {
            SWAP(i1,j1);
          }

          double value = stVKReducedInternalForces->cubicCoef(output,i1,j1,k1);

          if ((i == j) && (i == k)) // q_i^3
            value *= 3;
          else if ((i == j) || (i == k)) // q_i^2 * q_j
            value *= 2;

          quadraticCoef_[quadraticCoefPos(output,i,j,k)] = value;
        }
  }

  if (verbose)
    printf("\n");
  
  InitBuffers();
}


void StVKReducedStiffnessMatrix::Evaluate(double * q, double * Rq)
{
  // this is same as EvaluateSubset with start=0, end=quadraticSize

  /*
  int i,j,k;
  int output;

  // reset to free terms
  int index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      Rq[indexEntry] = freeCoef_[index];
      index++;
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // add linear terms
  index = 0;
  indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      for(j=0; j<r; j++)
      {
        Rq[indexEntry] += linearCoef_[index] * q[j];
        index++;
      }
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // add quadratic terms
  index = 0;
  indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      for(j=0; j<r; j++)
        for(k=j; k<r; k++)
        {
          Rq[indexEntry] += quadraticCoef_[index] * q[j] * q[k];
          index++;
        }
        indexEntry++;
    }
    indexEntry += output + 1;
  }

  // make symetric
  for(output=0; output<r; output++)
    for(i=0; i<output; i++)
      Rq[ELT(r,i,output)] = Rq[ELT(r,output,i)];
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

  // reset to free terms
  memcpy(buffer1,freeCoef_,sizeof(double)*quadraticSize);

  // add linear terms
  // multiply linearCoef_ and q
  // linearCoef_ is r x quadraticSize array
  cblas_dgemv(CblasColMajor, CblasTrans, 
        r, quadraticSize,
        1.0,
        linearCoef_, r,
        q, 1,
        1.0,
        buffer1, 1);

  // compute qiqj
  int index = 0;
  for(int output=0; output<r; output++)
    for(int i=output; i<r; i++)
    {
      qiqj[index] = q[output] * q[i];
      index++;
    }
 
  // update Rq
  // quadraticCoef_ is quadraticSize x quadraticSize matrix
  // each column gives quadratic coef for one matrix entry
  cblas_dgemv(CblasColMajor, CblasTrans, 
        quadraticSize, quadraticSize,
        1.0,
        quadraticCoef_, quadraticSize,
        qiqj, 1,
        1.0,
        buffer1, 1);

  // unpack into a symmetric matrix
  int i1=0,j1=0;
  for(int i=0; i< quadraticSize; i++)
  {
    Rq[ELT(r,i1,j1)] = buffer1[i];
    Rq[ELT(r,j1,i1)] = buffer1[i];
    j1++;
    if(j1 == r)
    {
      i1++;
      j1 = i1;
    }
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

void StVKReducedStiffnessMatrix::EvaluateSubset
  (double * q, int start, int end, double * Rq)
{
  /*
  int i,j,k;
  int output;

  // reset to free terms
  int index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      Rq[indexEntry] = freeCoef_[index];
      index++;
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // add linear terms
  index = 0;
  indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      for(j=0; j<r; j++)
      {
        Rq[indexEntry] += linearCoef_[index] * q[j];
        index++;
      }
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // add quadratic terms
  index = 0;
  indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      for(j=0; j<r; j++)
        for(k=j; k<r; k++)
        {
          Rq[indexEntry] += quadraticCoef_[index] * q[j] * q[k];
          index++;
        }
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // make symetric
  for(output=0; output<r; output++)
    for(i=0; i<output; i++)
      Rq[ELT(r,i,output)] = Rq[ELT(r,output,i)];
  */

  // reset to free terms
  memcpy(buffer1 + start,freeCoef_ + start,sizeof(double)*(end-start));

  // add linear terms
  // multiply linearCoef_ and q
  // linearCoef_ is r x quadraticSize array
  cblas_dgemv(CblasColMajor, CblasTrans, 
      r, end-start,
      1.0,
      linearCoef_ + r * start, r,
      q, 1,
      1.0,
      buffer1 + start, 1);

  // compute qiqj
  int index = 0;
  for(int output=0; output<r; output++)
    for(int i=output; i<r; i++)
    {
      qiqj[index] = q[output] * q[i];
      index++;
    }
 
  // update Rq
  // quadraticCoef_ is quadraticSize x quadraticSize matrix
  // each column gives quadratic coef for one matrix entry
  cblas_dgemv(CblasColMajor, CblasTrans, 
        quadraticSize, end-start,
        1.0,
        quadraticCoef_ + quadraticSize * start, quadraticSize,
        qiqj, 1,
        1.0,
        buffer1 + start, 1);

  // unpack into a symmetric matrix

  // first: set i1,j1
  int i1=0,j1=0;
  int iter=0;
  while (iter + r - i1 <= start)
  {
    iter += r - i1;
    i1++;
  }
  j1 = i1 + start - iter;
  // now (i1, j1) corresponds to i=start

  //if (matrixPosition(i1,j1) != start)
  //{
    //printf("Error.\n");
    //exit(1);
  //}

  for(int i=start; i< end; i++)
  {
    Rq[ELT(r,i1,j1)] = buffer1[i];
    Rq[ELT(r,j1,i1)] = buffer1[i];
    j1++;
    if(j1 == r)
    {
      i1++;
      j1 = i1;
    }
  }
}



void StVKReducedStiffnessMatrix::EvaluateLinear(double * q, double * Rq)
{
  int i;
  int output;

  // reset to free terms
  int index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    for(i=output; i<r; i++)
    {
      Rq[indexEntry] = freeCoef_[index];
      index++;
      indexEntry++;
    }
    indexEntry += output + 1;
  }

  // make symmetric
  indexEntry = 0;
  for(output=0; output<r; output++)
    for(i=0; i<output; i++)
      Rq[output * r + i] = Rq[i * r + output];

}

int StVKReducedStiffnessMatrix::Save(const char * filename)
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

  if ((int)(fwrite(freeCoef_,sizeof(double),r*(r+1)/2,fout)) < r*(r+1)/2)
    return 1;

  if ((int)(fwrite(linearCoef_,sizeof(double),r*(r+1)/2*linearSize,fout)) < r*(r+1)/2*linearSize)
    return 1;

  if ((int)(fwrite(quadraticCoef_,sizeof(double),r*(r+1)/2*quadraticSize,fout)) < r*(r+1)/2*quadraticSize)
    return 1;

  fclose(fout);

  return 0;
}


StVKReducedStiffnessMatrix::StVKReducedStiffnessMatrix(const char * filename)
{
  FILE * fin = fopen(filename,"rb");

  if (!fin)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  if ((int)(fread(&r,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  r2 = r * r;

  if ((int)(fread(&linearSize,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  if ((int)(fread(&quadraticSize,sizeof(int),1,fin)) < 1)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  freeCoef_ = (double*) malloc (sizeof(double) * r*(r+1)/2);

  if ((int)(fread(freeCoef_,sizeof(double),r*(r+1)/2,fin)) < r*(r+1)/2)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  linearCoef_ = (double*) malloc (sizeof(double) * r*(r+1)/2 * linearSize);

  if ((int)(fread(linearCoef_,sizeof(double),r*(r+1)/2*linearSize,fin)) < r*(r+1)/2*linearSize)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  quadraticCoef_ = (double*) malloc (sizeof(double) * r*(r+1)/2 * quadraticSize);

  if ((int)(fread(quadraticCoef_,sizeof(double),r*(r+1)/2*quadraticSize,fin)) < r*(r+1)/2*quadraticSize)
  {
    printf("Error: couldn't read from input stiffness matrix file.\n");
    return;
  }

  fclose(fin);

  InitBuffers();
}


void StVKReducedStiffnessMatrix::PrintMatrix()
{
  int i,j,k;
  int output;

  // free terms
  int index = 0;
  printf("{");
  for(output=0; output<r; output++)
  {
    printf("{");
    for(i=0; i<r; i++)
    {
      int i1 = output;
      int j1 = i;
      if (i1 > j1)
      {
        i1 = i;
        j1 = output;
      }
      printf("%.15f",freeCoef(i1,j1));
      if (i != r-1)
        printf(", ");
      index++;
      if (index % 5 == 4)
        printf("\n");
    }
    printf("}");
    if (output != r-1)
      printf(", ");
  }

  printf("} + ");

  // linear terms
  printf("{");
  index = 0;
  int indexEntry = 0;
  for(output=0; output<r; output++)
  {
    printf("{");
    for(i=0; i<r; i++)
    {
      int i1 = output;
      int j1 = i;
      if (i1 > j1)
      {
        i1 = i;
        j1 = output;
      }

      for(j=0; j<r; j++)
      {
        printf("%.15f * q%d",linearCoef(i1,j1,j),j);
        if (j != r - 1)
          printf(" + ");
        index++;

        if (index % 5 == 4)
          printf("\n");
      }

      if (i != r - 1)
        printf(", ");

      indexEntry++;
    }
    printf("}");
    if (output != r - 1)
      printf(",\n");
  }

  printf("}");

  // quadratic terms
  printf(" + {");
  index = 0;
  indexEntry = 0;
  for(output=0; output<r; output++)
  {
    printf("{");
    for(i=0; i<r; i++)
    {
      int i1 = output;
      int j1 = i;
      if (i1 > j1)
      {
        i1 = i;
        j1 = output;
      }

      for(j=0; j<r; j++)
        for(k=j; k<r; k++)
        {
          printf("%.15f * q%d * q%d",quadraticCoef(i1,j1,j,k),j,k);
          if (!((k == r - 1) && ( j == r-1)))
            printf(" + ");
          index++;

          if (index % 5 == 4)
            printf("\n");
        }

        if (i != r - 1)
          printf(", ");

        indexEntry++;
    }

    printf("}");
    if (output != r - 1)
      printf(",\n");
  }

  printf("}");
}

void StVKReducedStiffnessMatrix::UseSingleThread(int useSingleThread_)
{
  useSingleThread = useSingleThread_;
}

void StVKReducedStiffnessMatrix::InitBuffers()
{
  qiqj = (double*) malloc (sizeof(double) * quadraticSize);
  buffer1 = (double*) malloc (sizeof(double) * quadraticSize);
}

void StVKReducedStiffnessMatrix::FreeBuffers()
{
  free(qiqj);
  free(buffer1);
}

StVKReducedStiffnessMatrix * StVKReducedStiffnessMatrix::ShallowClone()
{
  StVKReducedStiffnessMatrix * output = new StVKReducedStiffnessMatrix(*this); // invoke default copy constructor
  output->shallowCopy = 1;
  output->InitBuffers();
  return output;
}

