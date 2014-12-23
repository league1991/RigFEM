/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "objMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC        *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder      *
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

#ifndef __BOUNDINGBOX_H__
#define __BOUNDINGBOX_H__

//  Bounding Box
//  Author: Jernej Barbic, CMU

#ifdef WIN32
  //#include <windows.h>
#endif

#include <vector>
#include "minivector.h"

class BoundingBox
{
public:

  BoundingBox(Vec3d bmin_g=Vec3d(0,0,0), Vec3d bmax_g=Vec3d(1,1,1)): bmin_(bmin_g), bmax_(bmax_g) { updateData();}
  template<class Triangle> BoundingBox(std::vector<Triangle> & tripool);

  // accessors
  Vec3d bmin() { return bmin_;}
  Vec3d bmax() { return bmax_;}

  Vec3d center() { return center_;}
  Vec3d halfSides() { return halfSides_;}

  double diameter() { return len(bmax_-bmin_); }

  // mutators
  void setbmin(Vec3d bmin_g) { bmin_=bmin_g; updateData();}
  void setbmin(double x, double y, double z) { bmin_=Vec3d(x,y,z); updateData();}
  void setbmax(Vec3d bmax_g) { bmax_=bmax_g; updateData();}
  void setbmax(double x, double y, double z) { bmax_=Vec3d(x,y,z); updateData();}

  void render();

  double distanceToPoint(Vec3d point);

  // sanity check bmin <= bmax
  void verifyBox();

  // expands from the center 
  // factor of 1.0 indicates no expansion
  void expand(double expansionFactor);
  void regularize(); // converts the box into one with all sides equal

  bool lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersection);

  void print();

protected:

  void updateData(); // updates center and half-sides
  Vec3d center_, halfSides_;
  Vec3d bmin_,bmax_;
};

#endif
