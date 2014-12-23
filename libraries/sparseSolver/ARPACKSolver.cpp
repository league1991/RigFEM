/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "ARPACKSolver.h"
#include "invMKSolver.h"
#include "lapack-headers.h"

//#define SPOOLES
#define PARDISO

#define DSAUPD DSAUPD_
#define DSEUPD DSEUPD_
#define INTEGER int

/* Note: 
  This code is now independent of ARPACK++ and calls ARPACK directly. It has been tested and works fine.

  Previously, it used ARPACK++, which, while convenient, required installing ARPACK++ (which required
  editing some of its headers files to bring it up to date wit today's C compilers), and
  only worked under Windows and Linux. Under Mac, it lead to strange errors inside ARPACK.
*/

// always keep this undefined
//#define USEARPACKPLUSPLUS

#ifdef USEARPACKPLUSPLUS
  #include "argsym.h"
#else
  extern "C"
  {
  void DSAUPD(INTEGER * IDO, char * BMAT, INTEGER * N, char * WHICH, INTEGER * NEV,
    double * TOL, double * RESID, INTEGER * NCV, double * V, INTEGER * LDV, INTEGER * IPARAM,
    INTEGER * IPNTR, double * WORKD, double * WORKL, INTEGER * LWORKL, INTEGER * INFO);

  // call DSAUPD 
  //c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
  //c       IPNTR, WORKD, WORKL, LWORKL, INFO )

  void DSEUPD(
    INTEGER * RVEC, char * HOWMNY, INTEGER * SELECT, double * D, double * Z, INTEGER * LDZ, double * SIGMA,
    char * BMAT, INTEGER * N, char * WHICH, INTEGER * NEV,
    double * TOL, double * RESID, INTEGER * NCV, double * V, INTEGER * LDV, INTEGER * IPARAM,
    INTEGER * IPNTR, double * WORKD, double * WORKL, INTEGER * LWORKL, INTEGER * INFO);

  //c  call DSEUPD  
  //c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
  //c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
}
#endif

#ifdef SPOOLES
  #include "SPOOLESSolver.h"
  #include "SPOOLESSolverMT.h"
#endif

#ifdef PARDISO
  #include "PardisoSolver.h"
#endif

int ARPACKSolver::SolveGenEigReg(SparseMatrix * KInput, SparseMatrix * MInput, int numEigenvalues, double * eigenvalues, double * eigenvectors, char * mode, int numLinearSolverThreads, int verbose)
{
  int np = KInput->Getn();
  if (np != MInput->Getn())
    return -1;

  if ((strcmp(mode, "SM") != 0) && (strcmp(mode, "LM") != 0))
    return -1;

  bool largestMagnitude = (strcmp(mode, "LM") == 0);
  SparseMatrix * K;
  SparseMatrix * M;
  if (largestMagnitude)
  {
    K = KInput;
    M = MInput;
  }
  else
  {
    // swap K and M to get the smallest eigenvalues
    K = MInput;
    M = KInput;
  }

  // solve Kx = lambda Mx with ARPACK
  // need multiplication with M^{-1} K and with M

  // create M^{-1} solver
  LinearSolver * invMSolver = NULL;
  #ifdef SPOOLES
    if (numLinearSolverThreads == 0)
    {
      if (verbose >= 1)
        printf("Linear solver: SPOOLES\n");
      invMSolver = new SPOOLESSolver(M, verbose);
    }
    else
    {
      if (verbose >= 1)
        printf("Linear solver: SPOOLES (%d threads).\n", numLinearSolverThreads);
      invMSolver = new SPOOLESSolverMT(M, numLinearSolverThreads, verbose);
    }
  #endif

  #ifdef PARDISO
    if (verbose >= 1)
      printf("Linear solver: PARDISO (%d threads).\n", (numLinearSolverThreads == 0) ? 1 : numLinearSolverThreads);
    int positiveDefinite = 0; 
    int directIterative = 0;
    PardisoSolver * pardisoSolver = new PardisoSolver(M, numLinearSolverThreads, positiveDefinite, directIterative, verbose);
    pardisoSolver->ComputeCholeskyDecomposition(M);
    invMSolver = pardisoSolver;
  #endif

  if (invMSolver == NULL)
  {
    printf("Error: no linear solver specified in SolveGenEigReg.\n");
    return 1;
  }

  // create M^{-1} K solver
  InvMKSolver invMKSolver(invMSolver, K);

  #ifdef USEARPACKPLUSPLUS
/*
    ARSymGenEig(int np, int nevp, ARFOP* objOPp,
              void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]), ARFB* objBp,
              void (ARFB::* MultBxp)(ARFLOAT[], ARFLOAT[]),
              char* whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
              int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
*/
    ARSymGenEig<double, InvMKSolver, SparseMatrix> solver
      (np, numEigenvalues, 
       &invMKSolver, &InvMKSolver::ComputeInvMK,
       M, (void (SparseMatrix::*)(double*,double*)) &SparseMatrix::MultiplyVector, "LM");

    int nconv = solver.EigenValVectors(eigenvectors, eigenvalues);
  #else
    const int maxIter = 10000;

    // call ARPACK
    INTEGER IDO = 0;
    char BMAT = 'G';
    INTEGER N = np;
    char WHICH[3] = "LM"; 
    INTEGER NEV = numEigenvalues;
    double TOL = 0.0;
    double * RESID = (double*) malloc (sizeof(double) * N);
    INTEGER NCV = 3 * NEV;
    if (NCV > N)
      NCV = N;
    double * V = (double*) malloc (sizeof(double) * N * NCV);
    INTEGER LDV = N;
    INTEGER IPARAM[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    IPARAM[0] = 1; // no user-provided shifts
    IPARAM[2] = maxIter; 
    IPARAM[3] = 1;
    //IPARAM[4] = 0;
    IPARAM[6] = 2; // the mode (regular generalized eigenproblem)
    INTEGER IPNTR[11]; 
    double * WORKD = (double*) malloc (sizeof(double) * 3 * N);
    INTEGER LWORKL = 2 * (NCV*NCV + 8*NCV); // at least NCV**2 + 8*NCV; 
    double * WORKL = (double*) malloc (sizeof(double) * LWORKL);
    INTEGER INFO = 0;

    double * buffer = (double*) malloc (sizeof(double) * N);
  
    do
    {
      //printf("Entering dsaupd...\n"); fflush(NULL);
      DSAUPD(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
      if (INFO != 0)
      {
        printf("Error: DSAUPD returned non-zero exit code %d .\n", (int)INFO);
        delete(invMSolver);
        return -1;
      }

      //printf("IDO=%d\n", (int)IDO); fflush(NULL);
  
      switch(IDO)
      {
        case -1:
          invMKSolver.ComputeInvMK(&WORKD[IPNTR[0]-1], &WORKD[IPNTR[1]-1]);
        break;

        case 1:
          invMKSolver.ComputeInvMK(&WORKD[IPNTR[0]-1], &WORKD[IPNTR[1]-1]);
        break;

        case 2:
          M->MultiplyVector(&WORKD[IPNTR[0]-1], &WORKD[IPNTR[1]-1]);
        break;

        case 3:
          printf("Error: case IDO=3 should have never happened.\n");
          delete(invMSolver);
          return -1;
        break;

        case 99:
        break;

        default:
          printf("Error: unknown case.\n");
          delete(invMSolver);
          return -1;
        break;
      };
    }
    while (IDO != 99);

    // obtain the eigenvalues and eigenvectors
    INTEGER RVEC = 1;
    char HOWMNY = 'A'; // all eigenvectors
    INTEGER * SELECT = (INTEGER*) malloc (sizeof(INTEGER) * NCV);
    double * D = eigenvalues;
    double * Z = eigenvectors;
    INTEGER LDZ = N;
    double * SIGMA = NULL;
    DSEUPD(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, SIGMA,
      &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
  
    free(SELECT); 
    free(RESID);
    free(V);
    free(WORKD);
    free(WORKL);
    free(buffer);

    int nconv = IPARAM[4];
  #endif

  // normalize results
  for(int i=0; i<nconv; i++)
    MInput->NormalizeVector(&eigenvectors[np*i]);

  //solver.FindEigenvectors();
  //nconv = solver.ConvergedEigenvalues();

  if (!largestMagnitude)
  {
    for(int i=0; i<nconv; i++)
      eigenvalues[i] = 1.0 / eigenvalues[i];

    // reverse eigenvalues (and appropriate eigenvectors)
    for(int i=0; i< numEigenvalues / 2; i++)
    {
      // swap i and numEigenvalues-1-i
      double temp = eigenvalues[i];
      eigenvalues[i] = eigenvalues[numEigenvalues-1-i];
      eigenvalues[numEigenvalues-1-i] = temp;

      for(int j=0; j<np; j++)
      {
        temp = eigenvectors[np*i+j];
        eigenvectors[np*i+j] = eigenvectors[np*(numEigenvalues-1-i) + j];
        eigenvectors[np*(numEigenvalues-1-i) + j] = temp;
      }
    }
  }

  if (verbose >= 1)
  {
    #ifdef USEARPACKPLUSPLUS
      cout << "Dimension of the system            : " << np              << endl;
      cout << "Number of 'requested' eigenvalues  : " << solver.GetNev()  << endl;
      cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
      cout << "Number of Arnoldi vectors generated: " << solver.GetNcv()  << endl;
      cout << "Number of iterations taken         : " << solver.GetIter() << endl;
    #else
      cout << "Dimension of the system            : " << np              << endl;
      cout << "Number of 'requested' eigenvalues  : " << NEV << endl;
      cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
      cout << "Number of Arnoldi vectors generated: " << NCV << endl;
      cout << "Number of iterations taken         : " << IPARAM[2] << endl;
      cout << endl;
    #endif

    // Printing eigenvalues.

    cout << "Eigenvalues:" << endl;
    //for (int i=numEigenvalues-nconv; i<numEigenvalues; i++) 
    for (int i=0; i<nconv; i++) 
    {
      cout << "  lambda[" << (i+1) << "]: " << eigenvalues[i] << endl;
    }
    cout << endl;

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    double * Ax      = new double[np];
    double * Bx      = new double[np];
    double * ResNorm = new double[nconv];

    double infinityNormK = KInput->GetInfinityNorm();
    double infinityNormM = MInput->GetInfinityNorm();
    double infNorm = infinityNormK;
    if (infinityNormM > infNorm)
      infNorm = infinityNormM;

    for (int i=0; i<nconv; i++) 
    {
      KInput->MultiplyVector(&eigenvectors[np*i],Ax);
      MInput->MultiplyVector(&eigenvectors[np*i],Bx);
      for(int j=0; j<np; j++)
        Ax[j] = Ax[j] - eigenvalues[i] * Bx[j];
      //cblas_daxpy(np, -eigenvalues[i], Bx, 1, Ax, 1);
      ResNorm[i] = 0;
      for(int j=0; j<np; j++)
        ResNorm[i] += Ax[j] * Ax[j];
      ResNorm[i] = sqrt(ResNorm[i]) / fabs(eigenvalues[i]) / infNorm;
      //ResNorm[i] = cblas_dnrm2(np, Ax, 1) / fabs(eigenvalues[i]);
    }

/*
    for (int i=0; i<nconv; i++) 
    {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*B*x(" << (i+1) << ")|| / |lambda| / max(||A||,||B||) = " << ResNorm[i] << endl;
    }
    cout << endl;
*/

    //printf("Cleaning up ARPACK workspace.\n");fflush(NULL);
    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;
  }

  delete(invMSolver);

  return nconv;
}

int ARPACKSolver::SolveGenEigShInv(SparseMatrix * K, SparseMatrix * M, int numEigenvalues, double * eigenvalues, double * eigenvectors, double sigma, int numLinearSolverThreads, int verbose)
{
  int np = K->Getn();
  if (np != M->Getn())
    return -1;

  // solve Kx = lambda Mx with ARPACK, shift-invert mode (mode number 3)
  // need multiplication with K^{-1} M and with M

  SparseMatrix * KsigmaM = K;
  if (sigma != 0)
  {
    // compute KsigmaM = K - sigma * M
    KsigmaM = new SparseMatrix(*K); 
    KsigmaM->BuildSubMatrixIndices(*M);
    KsigmaM->AddSubMatrix(-sigma, *M);
  }

  // create (K-sigma*M)^{-1} solver
  LinearSolver * invKSolver = NULL;

  #ifdef SPOOLES
    if (numLinearSolverThreads == 0)
    {
      if (verbose >= 1)
        printf("Linear solver: SPOOLES\n");
      invKSolver = new SPOOLESSolver(KsigmaM, verbose);
    }
    else
    {
      if (verbose >= 1)
        printf("Linear solver: SPOOLES (%d threads).\n", numLinearSolverThreads);
      invKSolver = new SPOOLESSolverMT(KsigmaM, numLinearSolverThreads, verbose);
    }
  #endif

  #ifdef PARDISO
    if (verbose >= 1)
      printf("Linear solver: PARDISO (%d threads).\n", (numLinearSolverThreads == 0) ? 1 : numLinearSolverThreads);
    int positiveDefinite = 0; 
    int directIterative = 0;
    PardisoSolver * pardisoSolver = new PardisoSolver(KsigmaM, numLinearSolverThreads, positiveDefinite, directIterative, verbose);
    pardisoSolver->ComputeCholeskyDecomposition(KsigmaM);
    invKSolver = pardisoSolver;
  #endif

  if (invKSolver == NULL)
  {
    printf("Error: no linear solver specified in SolveGenEigShInv.\n");
    return 1;
  }

  // create (K-sigma*M)^{-1}*M solver
  InvMKSolver invKMSolver(invKSolver, M);

  const int maxIter = 10000;

  // call ARPACK
  INTEGER IDO = 0;
  char BMAT = 'G';
  INTEGER N = np;
  char WHICH[3] = "LM"; 
  INTEGER NEV = numEigenvalues;
  double TOL = 0.0;
  double * RESID = (double*) malloc (sizeof(double) * N);
  INTEGER NCV = 3 * NEV;
  if (NCV > N)
    NCV = N;
  double * V = (double*) malloc (sizeof(double) * N * NCV);
  INTEGER LDV = N;
  INTEGER IPARAM[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  IPARAM[0] = 1; // no user-provided shifts
  IPARAM[2] = maxIter; 
  IPARAM[3] = 1;
  IPARAM[6] = 3; // the mode (shift-inverted generalized eigenproblem)
  INTEGER IPNTR[11]; 
  double * WORKD = (double*) malloc (sizeof(double) * 3 * N);
  INTEGER LWORKL = 2 * (NCV*NCV + 8*NCV); // at least NCV**2 + 8*NCV; 
  double * WORKL = (double*) malloc (sizeof(double) * LWORKL);
  INTEGER INFO = 0;

  double * buffer = (double*) malloc (sizeof(double) * N);
  
  if (verbose >= 1)
  {
    cout << "Dimension of the system            : " << np  << endl;
    cout << "Number of 'requested' eigenvalues  : " << NEV << endl;
    cout << "Entering the ARPACK eigenvalue routine..." << endl;
	fflush(NULL);
  }

  do
  {
    DSAUPD(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
    if (INFO != 0)
    {
      printf("Error: DSAUPD returned non-zero exit code %d .\n", (int)INFO);
      delete(invKSolver);
      return -1;
    }

    //printf("IDO=%d\n", (int)IDO); fflush(NULL);
    //c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
  
    switch(IDO)
    {
      case -1:
//c          IDO = -1: compute  Y = OP * X  where
//c                    IPNTR(1) is the pointer into WORKD for X,
//c                    IPNTR(2) is the pointer into WORKD for Y.
//c                    This is for the initialization phase to force the
//c                    starting vector into the range of OP.
        //printf("IDO = -1\n");
        invKMSolver.ComputeInvMK(&WORKD[IPNTR[0]-1], &WORKD[IPNTR[1]-1]);
      break;

      case 1:
//c          IDO =  1: compute  Y = (K - sigma*M)^-1 * Z
//c                    IPNTR(2) is the pointer into WORKD for Y,
//c                    IPNTR(3) is the pointer into WORKD for Z.
//c           (see dsdrv4.f example)
        //printf("IDO = 1\n");
        invKSolver->SolveLinearSystem(&WORKD[IPNTR[1]-1], &WORKD[IPNTR[2]-1]);
      break;

      case 2:
//c          IDO =  2: compute  Y = B * X  where
//c                    IPNTR(1) is the pointer into WORKD for X,
//c                    IPNTR(2) is the pointer into WORKD for Y.
        //printf("IDO = 2\n");
        M->MultiplyVector(&WORKD[IPNTR[0]-1], &WORKD[IPNTR[1]-1]);
      break;

      case 3:
        printf("Error: case IDO=3 should have never happened.\n");
        delete(invKSolver);
        return -1;
      break;

      case 99:
      break;

      default:
        printf("Error: unknown case.\n");
        delete(invKSolver);
        return -1;
      break;
    };
  }
  while (IDO != 99);

  // obtain the eigenvalues and eigenvectors
  INTEGER RVEC = 1;
  char HOWMNY = 'A'; // all eigenvectors
  INTEGER * SELECT = (INTEGER*) malloc (sizeof(INTEGER) * NCV);
  double * D = eigenvalues;
  double * Z = eigenvectors;
  INTEGER LDZ = N;
  double SIGMA = sigma;
  DSEUPD(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, &SIGMA,
    &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO);
  
  free(SELECT); 
  free(RESID);
  free(V);
  free(WORKD);
  free(WORKL);
  free(buffer);

  if (sigma != 0)
  {
    delete(KsigmaM);
  }

  int nconv = IPARAM[4];

  // normalize results
  for(int i=0; i<nconv; i++)
    M->NormalizeVector(&eigenvectors[np*i]);

  //solver.FindEigenvectors();
  //nconv = solver.ConvergedEigenvalues();

  if (verbose >= 1)
  {
    cout << "ARPACK solver is done." << endl;
    //cout << "Dimension of the system            : " << np              << endl;
    //cout << "Number of 'requested' eigenvalues  : " << NEV << endl;
    cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
    cout << "Number of Arnoldi vectors generated: " << NCV << endl;
    cout << "Number of iterations taken         : " << IPARAM[2] << endl;
    cout << endl;

    // Printing eigenvalues.

    cout << "Eigenvalues:" << endl;
    //for (int i=numEigenvalues-nconv; i<numEigenvalues; i++) 
    for (int i=0; i<nconv; i++) 
    {
      cout << "  lambda[" << (i+1) << "]: " << eigenvalues[i] << endl;
    }
    cout << endl;

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    double * Ax      = new double[np];
    double * Bx      = new double[np];
    double * ResNorm = new double[nconv];

    double infinityNormK = K->GetInfinityNorm();
    double infinityNormM = M->GetInfinityNorm();
    double infNorm = infinityNormK;
    if (infinityNormM > infNorm)
      infNorm = infinityNormM;

    for (int i=0; i<nconv; i++) 
    {
      K->MultiplyVector(&eigenvectors[np*i],Ax);
      M->MultiplyVector(&eigenvectors[np*i],Bx);
      for(int j=0; j<np; j++)
        Ax[j] = Ax[j] - eigenvalues[i] * Bx[j];
      //cblas_daxpy(np, -eigenvalues[i], Bx, 1, Ax, 1);
      ResNorm[i] = 0;
      for(int j=0; j<np; j++)
        ResNorm[i] += Ax[j] * Ax[j];
      ResNorm[i] = sqrt(ResNorm[i]) / fabs(eigenvalues[i]) / infNorm;
      //ResNorm[i] = cblas_dnrm2(np, Ax, 1) / fabs(eigenvalues[i]);
    }

/*
    for (int i=0; i<nconv; i++) 
    {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*B*x(" << (i+1) << ")|| / |lambda| / max(||A||,||B||) = " << ResNorm[i] << endl;
    }
    cout << endl;
*/

    //printf("Cleaning up ARPACK workspace.\n");fflush(NULL);
    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;
  }

  delete(invKSolver);
  return nconv;
}

/*
DSAUPD

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaupd
c
c\Description: 
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi 
c  Iteration.  For symmetric problems this reduces to a variant of the Lanczos 
c  method.  This method has been designed to compute approximations to a 
c  few eigenpairs of a linear operator OP that is real and symmetric 
c  with respect to a real positive semi-definite symmetric matrix B, 
c  i.e.
c                   
c       B*OP = (OP')*B.  
c
c  Another way to express this condition is 
c
c       < x,OPy > = < OPx,y >  where < z,w > = z'Bw  .
c  
c  In the standard eigenproblem B is the identity matrix.  
c  ( A' denotes transpose of A)
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dsaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x, A symmetric 
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dsaupd 
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dsaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dsaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          (If Mode = 2 see remark 5 below)
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * Z  and Z = B * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y,
c                    IPNTR(3) is the pointer into WORKD for Z.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c          After the initialization phase, when the routine is used in 
c          either the "shift-and-invert" mode or the Cayley transform
c          mode, the vector B * X is already available and does not 
c          need to be recomputed in forming OP*X.
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          Specify which of the Ritz values of OP to compute.
c
c          'LA' - compute the NEV largest (algebraic) eigenvalues.
c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
c          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
c          'BE' - compute NEV eigenvalues, half from each end of the
c                 spectrum.  When NEV is odd, compute one more from the
c                 high end than from the low end.
c           (see remark 1 below)
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector. 
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x. (See remark 4 below).
c
c  V       Double precision N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsaupd for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dsaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          Note: IPNTR(8:10) is only referenced by dseupd. See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd if RVEC = .TRUE. See Remarks.
c          Note: IPNTR(8:10) is only referenced by dseupd. See Remark 2.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine dseupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV .
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatable.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.
c
c
c\Remarks
c  1. The converged Ritz values are always returned in ascending 
c     algebraic order.  The computed Ritz values are approximate
c     eigenvalues of OP.  The selection of WHICH should be made
c     with this in mind when Mode = 3,4,5.  After convergence, 
c     approximate eigenvalues of the original problem may be obtained 
c     with the ARPACK subroutine dseupd. 
c
c  2. If the Ritz vectors corresponding to the converged Ritz values
c     are needed, the user must call dseupd immediately following completion
c     of dsaupd. This is new starting with version 2.1 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL'
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
c     linear systems should be solved with L and L' rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L'z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that NCV > NEV.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.   The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c
c  5. If IPARAM(7) = 2 then in the Reverse communication interface the user
c     must do the following. When IDO = 1, Y = OP * X is to be computed.
c     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
c     must overwrite X with A*X. Y is then the solution to the linear set
c     of equations B*Y = A*X.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) shifts in locations: 
c     1   WORKL(IPNTR(11))           
c     2   WORKL(IPNTR(11)+1)         
c                        .           
c                        .           
c                        .      
c     NP  WORKL(IPNTR(11)+NP-1). 
c
c     The eigenvalues of the current tridiagonal matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------------


Chao Yang
11/7/1997
*/
