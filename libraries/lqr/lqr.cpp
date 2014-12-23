/*************************************************************************
 *                                                                       *
 * "lqr" library , Copyright (C) 2009 MIT                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Jovan Popovic                                *
 * Funding: Singapore-MIT GAMBIT Game Lab                                *
 * Version 1.0                                                           *
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

#include <fstream>
#include "lqr.h"
#include "matrixIO.h"
#include "matrixExp.h"

using namespace std;

LQR::LQR(int stateSize_, int controlSize_, double dt_,
      double positionCost_, double controlCost_, double terminalCost_,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs, int verbose):
  stateSize(stateSize_), controlSize(controlSize_), dt(dt_)
{
  int T = Fs.size();
  printf("Building LQR controller (diagonal spherical costs)...\n");
  printf("T = %d, dt = %G, totalTime = %G\n"
         "stateSize = %d, controlSize = %d\n"
         "positionCost = %G, controlCost = %G, terminalCost = %G\n", 
         T, dt, T * dt, stateSize, controlSize, 
         positionCost_, controlCost_, terminalCost_); 

  double * positionCost = (double*) calloc (stateSize_ * stateSize_, sizeof(double));
  double * controlCost = (double*) calloc (controlSize_ * controlSize_, sizeof(double));
  double * terminalCost = (double*) calloc (stateSize_ * stateSize_, sizeof(double));

  for(int i=0; i<stateSize_; i++)
    positionCost[ELT(stateSize_, i, i)] = positionCost_;

  for(int i=0; i<controlSize_; i++)
    controlCost[ELT(controlSize_, i, i)] = controlCost_;

  for(int i=0; i<stateSize_; i++)
    terminalCost[ELT(stateSize_, i, i)] = terminalCost_;

  Init(stateSize_, controlSize_, dt_,
      positionCost, controlCost, terminalCost,
      Fs, Gs, verbose);

  free(terminalCost);
  free(controlCost);
  free(positionCost);
}

LQR::LQR(int stateSize_, int controlSize_, double dt_,
      double * positionCost, double * controlCost, double * terminalCost,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs, int verbose):
  stateSize(stateSize_), controlSize(controlSize_), dt(dt_)
{
  int T = Fs.size();
  printf("Building LQR controller (general dense cost matrices)...\n");
  printf("T = %d, dt = %G, totalTime = %G\n"
         "stateSize = %d, controlSize = %d\n"
         "positionCost[0,0] = %G, positionCost[%d,%d] = %G\n"
         "controlCost[0,0] = %G, controlCost[%d,%d] = %G\n"
         "terminalCost[0,0] = %G, terminalCost[%d,%d] = %G\n"
         "gain matrices footprint: %G Mb\n",
         T, dt, T * dt, stateSize, controlSize,
         positionCost[0], stateSize-1, stateSize-1, positionCost[stateSize*stateSize-1],
         controlCost[0], controlSize-1, controlSize-1, controlCost[controlSize*controlSize-1],
         terminalCost[0], stateSize-1, stateSize-1, terminalCost[stateSize*stateSize-1],
         1.0 * stateSize * controlSize * T * sizeof(double) / 1024 / 1024); 
  Init(stateSize_, controlSize_, dt_,
      positionCost, controlCost, terminalCost,
      Fs, Gs, verbose);
}

void LQR::Init(int stateSize_, int controlSize_, double dt_,
      double * positionCost, double * controlCost, double * terminalCost,
      const std::vector<double*> & Fs, const std::vector<double*> & Gs, int verbose)
{
  int T = Fs.size();

  Ks.clear();
  for(int i=0; i<=T; i++)
    Ks.push_back(Matrix<double>(controlSize, stateSize));
  Ks[T] = Matrix<double>(controlSize,stateSize,0.0);

  Matrix<double> Q(stateSize, stateSize, positionCost);
  Matrix<double> R(controlSize, controlSize, controlCost);

  // set initial Pk (terminal cost)
  Matrix<double> Pk(stateSize, stateSize, terminalCost);

  printf("Integrating Riccati equation... (dt=%G) ", dt);
  for (int i = T-1; i >= 0; i--) 
  {
    if (verbose)
    {
      printf(":");
      fflush(NULL);
    }

    // compute Phik, Gammak, Qk, Mk and Rk at step i    
    Matrix<double> Phik(stateSize, stateSize);
    Matrix<double> Gammak(stateSize, controlSize);
    Matrix<double> Qk(stateSize, stateSize);
    Matrix<double> Mk(stateSize, controlSize);
    Matrix<double> Rk(controlSize, controlSize);

    BuildRiccatiMatrices(dt, Fs[i], Gs[i],
      Q, R,
      Phik, Gammak,
      Qk, Mk, Rk);

    /*
      Phik.Save("A");
      Gammak.Save("B");
      Qk.Save("Q");
      Mk.Save("M");
      Rk.Save("R");
      exit(1);
    */

    // transpose matrices
    Matrix<double> PhiTk = Transpose(Phik);
    Matrix<double> GammaTk = Transpose(Gammak);
    Matrix<double> MTk = Transpose(Mk);

    // Riccati update
    if (verbose)
    {
      printf(".");
      fflush(NULL);
    }
    Matrix<double> A = Rk + GammaTk * Pk * Gammak;

    int rank;
    double rcond = 1e-12;
    Matrix<double> Ainv = LeastSquareSolve(A, Matrix<double>(controlSize, 1.0), rcond, &rank);
    if (rank < controlSize)
      printf("Warning: leastSquareSolve detected degenerate rank: %d (rcond=%G).\n", rank, rcond);

    // the Riccati recursive formula
    Matrix<double> Ck = Ainv * (MTk + GammaTk * Pk * Phik);
    Ks[i] = (-1.0) * Ck;
    Pk = Qk + (PhiTk * Pk * Phik) - (Mk + PhiTk * Pk * Gammak) * Ck;
  }

  P0 = new Matrix<double>(Pk);
  //P0->Print();

  if (!Nan_check())
  {
    printf("*** Error: LQR gain matrices contain nan. ***\n");
    printf("*********************************************\n");
    printf("*********************************************\n");
    printf("*********************************************\n");
    printf("*********************************************\n");
    //exit(1);
  }

  printf(" done\n");
}

LQR::LQR(char * filename)
{
  ifstream fin(filename,ios::binary);
  if (!fin)
    throw 1;

  fin.read((char*)&stateSize,sizeof(int));
  fin.read((char*)&controlSize,sizeof(int));
  fin.read((char*)&dt,sizeof(double));
  int T;
  fin.read((char*)&T,sizeof(int));
  for(int i=0; i<T; i++)
  {
    Ks.push_back(Matrix<double>(controlSize, stateSize, 0.0));
    fin.read((char*)(Ks[i].GetData()),stateSize * controlSize * sizeof(double));
  }
  P0 = new Matrix<double>(stateSize, stateSize, 0.0);
  fin.read((char*)(P0->GetData()),stateSize * stateSize * sizeof(double));

  fin.close();

  printf("Loaded LQR model from %s.\n", filename);
  printf("stateSize=%d controlSize=%d dt=%G T=%d\n", 
   stateSize, controlSize, dt, T);
}

void LQR::Save(char * filename)
{
  ofstream fout(filename,ios::binary);
  fout.write((char*)&stateSize,sizeof(int));
  fout.write((char*)&controlSize,sizeof(int));
  fout.write((char*)&dt,sizeof(double));
  int T = Ks.size();
  fout.write((char*)&T,sizeof(int));
  for(int i=0; i<T; i++)
    fout.write((char*)(Ks[i].GetData()),stateSize * controlSize * sizeof(double));
  fout.write((char*)(P0->GetData()),stateSize * stateSize * sizeof(double));

  fout.close();
}

LQR::~LQR()
{
  delete(P0);
}

void LQR::ComputeControl(int timestepIndex, double * dx_, double * du_)
{
  // multiply du = C * dx
  // Ks[timestepIndex] is controlSize x stateSize
  Matrix<double> dx(stateSize, 1, dx_, false, false);
  Matrix<double> du = Ks[timestepIndex] * dx;
  memcpy(du_, du.GetData(), sizeof(double) * controlSize);
}

Matrix<double> & LQR::GetGainMatrix(int timestepIndex)
{
  return Ks[timestepIndex];
}

bool LQR::my_isnan(double x)
{
  return (x != x);
}

bool LQR::Nan_check()
{
  int len = stateSize*controlSize;
  for(unsigned int i=0; i<Ks.size(); i++)
  {
    double * K = Ks[i].GetData();
    for(int j=0; j<len; j++)
    {
      if (my_isnan(K[j]))
        return false;
    }
  }

  return true;
}

double LQR::ComputeInitialCost(double * dx)
{
  Matrix<double> dxM(stateSize, 1, dx, false, false);
  Matrix<double> costM = 0.5 * Transpose(dxM) * (*P0 * dxM);
  return costM(0,0);
}

// builds the matrices necessary for the Riccati recursive formula
// some of these matrices involve integrals of matrix exponentials, which are approximated using Simpson's integration
void LQR::BuildRiccatiMatrices(double dt, double * F, double * G,
    Matrix<double> & Q, Matrix<double> & R,
    Matrix<double> & Phik, Matrix<double> & Gammak,
    Matrix<double> & Qk, Matrix<double> & Mk, Matrix<double> & Rk)
{
  const int nSimpson = 10; // the number of subintervals (plus one) for Simpson integration
  const int ideg = 6; // the accuracy degree for the MatrixExp routine

  Matrix<double> FM(stateSize, stateSize, F);
  Matrix<double> FInvM = Inverse(FM);
  Matrix<double> GM(stateSize, controlSize, G);

  // compute Phis and Gammas arrays, which will hold Phi and Gamma values at the nodes of the subintervals
  double * Phis = (double*) calloc (stateSize * stateSize * (nSimpson+1), sizeof(double));
  double * Gammas = (double*) calloc (stateSize * controlSize * (nSimpson+1), sizeof(double));
  #define PHI(i) (&Phis[(i) * stateSize * stateSize])
  #define GAMMA(i) (&Gammas[(i) * stateSize * controlSize])

  // compute Phis
  // Phi(i) = exp(F * i / nSimpson * dt)
  for(int i=0; i<stateSize; i++)
    PHI(0)[ELT(stateSize,i,i)] = 1.0;
  for(int i=1; i<=nSimpson; i++)
    MatrixExp(stateSize, FM.GetData(), 1.0 * i / nSimpson * dt, PHI(i), ideg);

  // compute Gammas
  // Gamma(i) = F^{-1} * (exp(F * i / nSimpson * dt) - I) * G
  // note that F^{-1} and (exp(F * i / nSimpson * dt) - I) commute (can reverse multiplication order) as they are both polynomials in F and F^{-1} (by a theorem in algebra)
  for(int i=0; i<=nSimpson; i++)
  {
    Matrix<double> PhiM(stateSize, stateSize, PHI(i), false, false);
    Matrix<double> GammaM(stateSize, controlSize, GAMMA(i), false, false);
    GammaM = (PhiM - Matrix<double>(stateSize,1.0)) * FInvM * GM;
  }

  // build Simpson weights
  double * simpsonWeights = (double*) malloc (sizeof(double) * (nSimpson+1));
  for(int i=0; i<=nSimpson; i++)
  {
    if ((i==0) || (i==nSimpson))
      simpsonWeights[i] = 1.0;
    else
    {
      if (i % 2 == 0)
        simpsonWeights[i] = 2.0;
      else
        simpsonWeights[i] = 4.0;
    }
  }  
  for(int i=0; i<=nSimpson; i++)
    simpsonWeights[i] *= 1.0 * dt / nSimpson / 3;

  // Phik = exp(F * dt)
  // Gammak = F^{-1} * (exp(F * dt) - I) * G
  Phik = Matrix<double> (stateSize,stateSize, PHI(nSimpson));
  Gammak = Matrix<double> (stateSize, controlSize, GAMMA(nSimpson));

  // Qk = integral_{tau=0}^{tau=dt} Phi(tau)^T * Q(t_i + tau) * Phi(tau) dtau
  // Mk = integral_{tau=0}^{tau=dt} Phi(tau)^T * Q(t_i + tau) * Gamma(tau) dtau
  // Rk = integral_{tau=0}^{tau=dt} ( Gamma(tau)^T * Q(t_i + tau) * Gamma(tau) + R(t_i + tau) ) dtau
  Qk = Matrix<double>(stateSize,stateSize,0.0);
  Mk = Matrix<double>(stateSize,controlSize,0.0);
  Rk = Matrix<double>(controlSize,controlSize,0.0);
  for(int i=0; i<=nSimpson; i++)
  {
    Matrix<double> PhiM(stateSize, stateSize, PHI(i), false, false);
    Matrix<double> GammaM(stateSize, controlSize, GAMMA(i), false, false);
    Qk += simpsonWeights[i] * (Transpose(PhiM) * Q * PhiM); 
    Mk += simpsonWeights[i] * (PhiM * Q * GammaM); 
    Rk += simpsonWeights[i] * (Transpose(GammaM) * Q * GammaM + R);
  }

  free(simpsonWeights);
  free(Phis);
  free(Gammas);
}

