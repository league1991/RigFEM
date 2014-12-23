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
  The driver for the gradient of polar decomposition of a general 3x3 matrix
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

/*

  You can compile this driver, say, by:
  g++ -o polarDecompositionGradientDriver polarDecompositionGradientDriver.cpp polarDecompositionGradient.cpp polarDecomposition.cpp vec3d.cpp mat3d.cpp eig3.cpp

*/

#include "polarDecomposition.h"
#include "polarDecompositionGradient.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintVector3(double * vec, char * s)
{
  if (s != NULL)
    printf("%s\n", s);

  printf("%G %G %G\n", vec[0], vec[1], vec[2]);
}

void PrintMatrix3x3(double * mat, char * s)
{
  if (s != NULL)
    printf("%s\n", s);

  printf("%G %G %G\n", mat[0], mat[1], mat[2]);
  printf("%G %G %G\n", mat[3], mat[4], mat[5]);
  printf("%G %G %G\n", mat[6], mat[7], mat[8]);
}

int main()
{
  // all matrices are 3x3, stored row-major
  // random matrix for M (for testing)
  double M[9] = {-0.2807952341472597, 0.3115011511885256, -0.0485903249146864, 
                 -0.45015464324442633, -0.8543680801174494, 0.19492773318804713, 
                  0.36840026899411693, 0.8634092118662448, 0.7528384426032098};     

  double Q[9]; // output
  double S[9]; // output
  double tol = 1e-12;
  PolarDecomposition::Compute(M, Q, S, tol); 

  PrintMatrix3x3(Q, "Q");
  PrintMatrix3x3(S, "S");

  // random matrix for MDot (for testing)
  double MDot[9] = {0.5118439624450675, 0.24235217728059955, -0.7450654689862747, 
                    -0.4345942647112335, 1.6156424648429977, -0.6830032318171733, 
		    -0.44285129704233783, -0.4366569686333789, -1.0401543827923478};

  // random matrix for MDotDot (for testing)
  double MDotDot[9] = {-0.12353971625289528, -0.14236178366865948, 0.7566400471581742, 
                       1.8378313423900419, 0.1496418600146152, 0.351961331124921, 
		       0.42477528592658276, 0.6794947974490568, 1.3110637123590698};

  double omega[3]; // output
  double QDot[9]; // output
  double SDot[9]; // output
  double omegaDot[3]; // output
  double QDotDot[9]; // output
  PolarDecompositionGradient::Compute(M, Q, S, MDot, omega, QDot, SDot, MDotDot, omegaDot, QDotDot);

  PrintVector3(omega, "omega");
  PrintMatrix3x3(QDot, "QDot");
  PrintMatrix3x3(SDot, "SDot");
  PrintVector3(omegaDot, "omegaDot");
  PrintMatrix3x3(QDotDot, "QDotDot");

  return 0;
}

