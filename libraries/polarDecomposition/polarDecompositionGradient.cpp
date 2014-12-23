/*

* Copyright (c) 2011, Jernej Barbic, Yili Zhao, University of Southern California
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of University of Southern California, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY JERNEJ BARBIC, YILI ZHAO AND UNIVERSITY OF SOUTHERN CALIFORNIA 
* ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
* IN NO EVENT SHALL JERNEJ BARBIC, YILI ZHAO OR UNIVERSITY OF SOUTHERN CALIFORNIA BE 
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
  Gradient of polar decomposition of a general 3x3 matrix
  Version 2.0

  This code was used in the following publication:

  Jernej Barbic, Yili Zhao:
  Real-time Large-deformation Substructuring, ACM Transactions on Graphics (SIGGRAPH 2011), 30(4), Aug 2011

  @article{Barbic:2011:RLS,
    author =  {Jernej Barbi\v{c} and Yili Zhao},
    journal = {ACM Trans. on Graphics (SIGGRAPH 2011)},
    number =  "4",
    title =   "Real-time Large-deformation Substructuring",
    volume =  "30",
    year =    "2011",
    pages =   "91:1--91:7",
  }

  Authors of this code: Jernej Barbic, Yili Zhao
*/

#include "polarDecompositionGradient.h"
#include "matrixMultiplyMacros.h"
#include "mat3d.h"

void PolarDecompositionGradient::Compute(const double * M, const double * Q, const double * S, const double * MDot, double * omega, double * QDot, double * SDot, const double * MDotDot, double * omegaDot, double * QDotDot)
{
  // compute omega = G^{-1} (2 * skew(Q^T * MDot)), where G = (tr(S)I - S) * Q^T
  // (see Barbic and Zhao, SIGGRAPH 2011)

  // first, construct G, and invert it

  // tempMatrix = tr(S)I - S
  double tempMatrix[9];
  for(int i=0; i<9; i++)
    tempMatrix[i] = -S[i];
  double trace = S[0] + S[4] + S[8];
  tempMatrix[0] += trace;
  tempMatrix[4] += trace;
  tempMatrix[8] += trace;

  double G[9]; // G = (tr(S)I - S) * Q^T
  MATRIX_MULTIPLY3X3ABT(tempMatrix, Q, G);
  Mat3d GM(G);
  Mat3d GInvM = inv(GM);
  double GInv[9];
  GInvM.convertToArray(GInv);

  // omega = GInv * (2 * skew(R^T * Mdot))
  MATRIX_MULTIPLY3X3ATB(Q, MDot, tempMatrix);
  double rhs[3];
  SKEW_PART(tempMatrix, rhs);
  VECTOR_SCALE3(rhs, 2.0);
  MATRIX_VECTOR_MULTIPLY3X3(GInv, rhs, omega);

  // compute QDot = tilde(omega) * Q
  double omegaTilde[9];
  SKEW_MATRIX(omega, omegaTilde);
  //double QDot[9];
  MATRIX_MULTIPLY3X3(omegaTilde, Q, QDot);

  // compute SDot = Q^T * (MDot - QDot * S)
  // tempMatrix = MDot - QDot * S
  MATRIX_MULTIPLY3X3(QDot, S, tempMatrix);
  for(int i=0; i<9; i++)
    tempMatrix[i] = MDot[i] - tempMatrix[i];
  // SDot = Q^T * tempMatrix
  MATRIX_MULTIPLY3X3ATB(Q, tempMatrix, SDot); 

  if ((MDotDot != NULL) && (omegaDot != NULL))
  {
    // compute omegaDot = GInv * ( 2 skew(Q^T (ADotDot - omegaTilde * ADot)) - (tr(SDot) I - SDot) * Q^T * omega )
    // (see Barbic and Zhao, SIGGRAPH 2011)
    
    // tempMatrix = MDotDot - omegaTilde * MDot
    MATRIX_MULTIPLY3X3(omegaTilde, MDot, tempMatrix);
    for(int i=0; i<9; i++)
      tempMatrix[i] = MDotDot[i] - tempMatrix[i];

    double tempMatrix2[9];
    // tempVector = 2 * skew(Q^T * tempMatrix)
    MATRIX_MULTIPLY3X3ATB(Q, tempMatrix, tempMatrix2);

    double tempVector[3];
    SKEW_PART(tempMatrix2, tempVector);
    VECTOR_SCALE3(tempVector, 2.0);

    // tempMatrix = tr(SDot)I - SDot
    for(int i=0; i<9; i++)
      tempMatrix[i] = -SDot[i];
    double trace = SDot[0] + SDot[4] + SDot[8];
    tempMatrix[0] += trace;
    tempMatrix[4] += trace;
    tempMatrix[8] += trace;

    // tempVector2 = (tempMatrix * Q^T) * omega
    double tempVector2[3];
    MATRIX_MULTIPLY3X3ABT(tempMatrix, Q, tempMatrix2);
    MATRIX_VECTOR_MULTIPLY3X3(tempMatrix2, omega, tempVector2);

    // tempVector -= tempVector2
    VECTOR_SUBTRACTEQUAL3(tempVector, tempVector2);

    // tempVector2 = GInv * tempVector
    MATRIX_VECTOR_MULTIPLY3X3(GInv, tempVector, omegaDot);

    if (QDotDot != NULL)
    {
      double tempMatrix[9];
      SKEW_MATRIX(omegaDot, tempMatrix);
      MATRIX_MULTIPLY3X3(omegaTilde, omegaTilde, tempMatrix2);
      for(int i=0;i<9;i++)
	tempMatrix[i] += tempMatrix2[i];
      MATRIX_MULTIPLY3X3(tempMatrix, Q, QDotDot);
    }
  }
}

