/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "quaternion" library , Copyright (C) 2007 CMU                         *
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

#include "quaternion.h"

/*
  See David Baraff's SIGGRAPH course notes for the description of the method below:
  "An Introduction to Physically Based Modeling:
   Rigid Body Simulation I: Unconstrained Rigid Body Dynamics"
  http://www.cs.cmu.edu/~baraff/pbm/pbm.html
*/

//Transforms the given matrix (assumed orthogonal, specified 
//in row-major order) into one of the two corresponding quaternions.
template <typename real>
Quaternion<real> Quaternion<real>::Matrix2Quaternion(real * R)
{
/* 
   Order of matrix elements is row-major:

   (0,0) 0  (0,1) 1  (0,2) 2
   (1,0) 3  (1,1) 4  (1,2) 5
   (2,0) 6  (2,1) 7  (2,2) 8
*/

  Quaternion<real> q;
  real trace, u;
  trace = R[0] + R[4] + R[8];

  if(trace >= 0)
  {
    u = (real)sqrt(trace + 1);
    q.s = (real)0.5 * u;
    u = (real)0.5 / u;
    q.x = (R[7] - R[5]) * u;
    q.y = (R[2] - R[6]) * u;
    q.z = (R[3] - R[1]) * u;
  }
  else
  {
    int i = 0;
    if(R[4] > R[0])
      i = 1;

    if(R[8] > R[3*i+i])
      i = 2;

    switch (i)
    {
      case 0:
        u = (real)sqrt((R[0] - (R[4] + R[8])) + 1);
        q.x = 0.5f * u;
        u = 0.5f / u;
        q.y = (R[3] + R[1]) * u;
        q.z = (R[2] + R[6]) * u;
        q.s = (R[7] - R[5]) * u;
      break;

      case 1:
        u = (real)sqrt((R[4] - (R[8] + R[0])) + 1);
        q.y = 0.5f * u;
        u = 0.5f / u;
        q.z = (R[7] + R[5]) * u;
        q.x = (R[3] + R[1]) * u;
        q.s = (R[2] - R[6]) * u;
      break;

      case 2:
        u = (real)sqrt((R[8] - (R[0] + R[4])) + 1);
        q.z = 0.5f * u;

        u = 0.5f / u;
        q.x = (R[2] + R[6]) * u;
        q.y = (R[7] + R[5]) * u;
        q.s = (R[3] - R[1]) * u;
      break;
    }
  }

  return q;
}

template Quaternion<double> Quaternion<double>::Matrix2Quaternion(double * R);
template Quaternion<float> Quaternion<float>::Matrix2Quaternion(float * R);

