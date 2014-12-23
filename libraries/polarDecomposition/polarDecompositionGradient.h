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

#ifndef _POLARDECOMPOSITIONGRADIENT_H_
#define _POLARDECOMPOSITIONGRADIENT_H_

#include <stdlib.h>

class PolarDecompositionGradient
{
public:
  // computes the first and second *derivative* of the rotation matrix in polar decomposition
  // also computes the first derivative of the symmetric matrix in polar decomposition
  // for more details, please see the citation above

  // Polar decomposition is:
  // M = Q * S
  // all matrices are 3x3, stored row major as a 9-vector of doubles

  // Gradient computation instructions:
  // inputs: M, Q, S, MDot, MDotDot (all 3x3 matrices)
  // the routine assumes that the polar decomposition has been already performed (e.g., using polarDecomposition.h)
  //   M: the input matrix (should be non-singular)
  //   Q: the polar decomposition rotation (obtained, say, via polarDecomposition.h)
  //   S: the polar decomposition symmetric matrix (obtained, say, via polarDecomposition.h)
  //   MDot: the derivative of M
  //   MDotDot: second derivative of M (only needed for omegaDot and QDotQDot, can be ommitted otherwise)
  // outputs: omega, QDot, SDot, omegaDot, QDotDot (all 3x3 matrices, except omega and omegaDot which are 3-vectors)
  //   omega: the rotational velocity of Q 
  //   QDot: derivative of Q  
  //   SDot: derivative of S
  //   omegaDot: derivative of omega (only computed when MDotDot is provided, not referenced otherwise)   
  //   QDotDot: second derivative of Q (only computed when MDotDot is provided and QDotDot is not NULL, not referenced otherwise)
  // See also the included example driver.
  static void Compute(const double * M, const double * Q, const double * S, const double * MDot, double * omega, double * QDot, double * SDot, const double * MDotDot=NULL, double * omegaDot=NULL, double * QDotDot = NULL);

protected:

};

#endif

