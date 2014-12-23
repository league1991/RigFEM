#include <stdio.h>
#include <math.h>
#include "polarDecomposition.h"
#include "mat3d.h"

/*
  See polarDecomposition.h for license information.
*/

// compute the one-norm of a 3x3 matrix (row-major)
double PolarDecomposition::oneNorm(const double * A)
{
  double norm = 0.0;
  for (int i=0; i<3; i++) 
  {
    double columnAbsSum = fabs(A[i + 0]) + fabs(A[i + 3]) + fabs(A[i + 6]);
    if (columnAbsSum > norm) 
      norm = columnAbsSum;
  }
  return norm;
}

// compute the inf-norm of a 3x3 matrix (row-major)
double PolarDecomposition::infNorm(const double * A)
{
  double norm = 0.0;
  for (int i=0; i<3; i++) 
  {
    double rowSum = fabs(A[3 * i + 0]) + fabs(A[3 * i + 1]) + fabs(A[3 * i + 2]);
    if (rowSum > norm) 
      norm = rowSum;
  }
  return norm;
}

// Input: M (3x3 mtx)
// Output: Q (3x3 rotation mtx), S (3x3 symmetric mtx)
double PolarDecomposition::Compute(const double * M, double * Q, double * S, double tolerance, int forceRotation)
{
  double Mk[9];
  double Ek[9];
  double det, M_oneNorm, M_infNorm, E_oneNorm;
  int useSVD = 0;

  // Mk = M^T
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      Mk[3 * i + j] = M[3 * j + i];

  M_oneNorm = oneNorm(Mk); 
  M_infNorm = infNorm(Mk);

  do 
  {
    double MadjTk[9];
 
    // row 2 x row 3
    crossProduct(&(Mk[3]), &(Mk[6]), &(MadjTk[0])); 
    // row 3 x row 1
    crossProduct(&(Mk[6]), &(Mk[0]), &(MadjTk[3]));
    // row 1 x row 2
    crossProduct(&(Mk[0]), &(Mk[3]), &(MadjTk[6]));

    det = Mk[0] * MadjTk[0] + Mk[1] * MadjTk[1] + Mk[2] * MadjTk[2];

    if ((det <= 1e-6) && forceRotation)
    {
      useSVD = 1;
      break;
    }

    if (det == 0.0) 
    {
      printf("Warning (polarDecomposition) : zero determinant encountered.\n");
      break;
    }

    double MadjT_one = oneNorm(MadjTk); 
    double MadjT_inf = infNorm(MadjTk);

    double gamma = sqrt(sqrt((MadjT_one * MadjT_inf) / (M_oneNorm * M_infNorm * det * det)));
    double g1 = gamma * 0.5;
    double g2 = 0.5 / (gamma * det);

    for(int i=0; i<9; i++)
    {
      Ek[i] = Mk[i];
      Mk[i] = g1 * Mk[i] + g2 * MadjTk[i];
      Ek[i] -= Mk[i];
    }

    E_oneNorm = oneNorm(Ek);
    M_oneNorm = oneNorm(Mk);  
    M_infNorm = infNorm(Mk);
  }
  while ( E_oneNorm > M_oneNorm * tolerance );

  if (useSVD)
  {
    // use the SVD algorithm to compute Q
    Mat3d Mm(M);
    double modifiedSVD_singularValue_eps = tolerance;
    Mat3d Um, Vm;
    Vec3d Lambda;
    int modifiedSVD = 1;
    int code = SVD(Mm, Um, Lambda, Vm, modifiedSVD_singularValue_eps, modifiedSVD);
    code = code + 1; // just to ignore compiler warning

    Mat3d Qm = Um * trans(Vm);
    Qm.convertToArray(Q);
  }
  else
  {
    // Q = Mk^T 
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        Q[3*i+j] = Mk[3*j+i];
  }

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      S[3*i+j] = 0;
      for(int k=0; k<3; k++)
        S[3*i+j] += Mk[3*i+k] * M[3*k+j];
    }
    
  // S must be symmetric; enforce the symmetry
  S[1] = S[3] = 0.5 * (S[1] + S[3]);
  S[2] = S[6] = 0.5 * (S[2] + S[6]);
  S[5] = S[7] = 0.5 * (S[5] + S[7]);

  return det;
}

