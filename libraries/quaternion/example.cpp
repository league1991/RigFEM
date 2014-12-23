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
#include "stdio.h"

#define PI 3.141592653589793238462643

// This is example code that demonstrates how to use the Quaternion class.
// Example author: Jernej Barbic

int main()
{
  // make a quaternion
  Quaternion<float> q(2.0, 1.5, 1.0, -0.5);
  // print it out to screen
  q.Print();

  // compute the norm
  float norm = q.Norm();
  printf("Quaternion norm is: %G\n", norm);

  // normalize it
  q.Normalize();
  printf("Normalized quaternion is:\n");
  q.Print(); 
  printf("===================\n");

  // OK, done with q, let's make some more quaternions

  // make a quaternion corresponding to a rotation of 90 degrees around x-axis
  float axis1[3] = {1, 0, 0};
  float angle1 = PI * 0.5; // 90 deg in radians
  Quaternion<float> q1(angle1, axis1);
  printf("First quaternion:\n");
  q1.Print();

  // make a quaternion corresponding to a rotation of 90 degrees around z-axis
  float axis2[3] = {0, 0, 1};
  float angle2 = PI * 0.5; // 90 deg in radians
  Quaternion<float> q2(angle2, axis2);
  printf("Second quaternion:\n");
  q2.Print();

  // compute a quaternion corresponding to a rotation composition:
  //   first, around x-axis, then around z-axis
  Quaternion<float> q3 = q2 * q1;
  printf("Quaternion corresponding to matrix composition is:\n");
  q3.Print();

  // find out the angle and axis of rotation corresponding to q3
  float angle3;
  float axis3[3];
  q3.GetRotation(&angle3, axis3);
  printf("Angle of rotation: %G deg. Axis: [%G %G %G].\n", angle3 * 180 / PI, axis3[0], axis3[1], axis3[2]);

  // covert q3 to the corresponding 3x3 rotation (given in row-major order)
  float R[9];
  q3.Quaternion2Matrix(R);

  // rotate a vector by R (multiply R by vector)
  float vector[3] = {0,0,1};
  float rotatedVector[3];
  for(int i=0; i<3; i++)
    rotatedVector[i] = R[3*i+0] * vector[0] + R[3*i+1] * vector[1] + R[3*i+2] * vector[2];
  printf("Original vector is: [%G %G %G].\n", vector[0], vector[1], vector[2]);
  printf("Rotated vector is: [%G %G %G].\n", rotatedVector[0], rotatedVector[1], rotatedVector[2]);

  // compute the inverse rotation
  //   for unit quaternions, inverse is simply the conjugate value
  Quaternion<float> q4 = q3.conj();
  float RInv[9];
  q4.Quaternion2Matrix(RInv);

  // sanity check: S = R * RInv must be the 3x3 identity matrix
  float S[9];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      S[3*i+j] = 0;
      for(int k=0; k<3; k++)
        S[3*i+j] += R[3*i+k] * RInv[3*k+j];
    }

  printf("The following should be the identity matrix (modulo numerical round-off):\n");
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
      printf("%G ", S[3*i+j]);
    printf("\n");
  }

  // convert S to quaternion
  Quaternion<float> q5 = Quaternion<float>::Matrix2Quaternion(S);
  // This will be either -1 or 1
  printf("Quaternion corresponding to the above unit matrix:\n");
  q5.Print();
    
  return 0;
}

