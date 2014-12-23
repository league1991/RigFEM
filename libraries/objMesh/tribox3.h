/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "objMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC        *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Somya Sharma, Jernej Barbic                             *
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

/*
  Tests if a 3D triangle overlaps with a 3D box.

  INPUT: center of box, box half sizes (in each of the three dimensions), and the three triangle vertices v0, v1, v2.
  OUTPUT: whether triangle intersects the box or not
  Note: all entries in boxHalfSize must be >= 0.
  Note: lower-left-front corner is boxCenter - boxHalfSize, upper-right-back corner is boxCenter + boxHalfSiz
e.
*/

#ifndef _TRIBOX3_H_
#define _TRIBOX3_H_

bool triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

#endif

