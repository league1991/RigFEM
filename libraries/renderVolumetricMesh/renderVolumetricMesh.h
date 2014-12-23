/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "renderVolumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT,    *
 *                                                          2014 USC     *
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

#ifndef _RENDER_VOLUMETRIC_MESH_
#define _RENDER_VOLUMETRIC_MESH_

#include <set>
#include "printBitmap.h"
#include "volumetricMesh.h"

class RenderVolumetricMesh
{
public:
  RenderVolumetricMesh();

  void Render(VolumetricMesh * volumetricMesh, int wireframe=0, double * u = NULL);
  void RenderWireframe(VolumetricMesh * volumetricMesh);
  void RenderSolidAndWireframe(VolumetricMesh * volumetricMesh);

  void RenderVertices(VolumetricMesh * volumetricMesh);
  void RenderVertices(VolumetricMesh * volumetricMesh, int * vertices, int numVertices, bool oneIndexed=true);
  void RenderVertices(VolumetricMesh * volumetricMesh, std::set<int> * vertices, bool oneIndexed=true);
  void SelectRenderVertices(VolumetricMesh * volumetricMesh);

  void RenderVertexLabels(VolumetricMesh * volumetricMesh);
  void RenderVertexLabels(VolumetricMesh * volumetricMesh, int start, int end);

  void DrawUnselectedPoints(VolumetricMesh * volumetricMesh, int * selectionArray);
  void DrawSelectedPoints(VolumetricMesh * volumetricMesh, int * selectedVertices, int numSelectedVertices);

  void RenderDeformation(VolumetricMesh * volumetricMesh, double * u);
  void RenderWireframeDeformation(VolumetricMesh * volumetricMesh, double * u);
  void RenderSolidAndWireframeDeformation(VolumetricMesh * volumetricMesh, double * u);
  void RenderVertexDeformed(VolumetricMesh * volumetricMesh, int ver, double * u);  

  // controls how different material groups are rendered
  void SetFlatRenderingMode(); // all colored white
  void SetGradedRenderingMode(VolumetricMesh * volumetricMesh); // colored graded, according to numerical values of E, nu, density
  void SetDiscreteRenderingMode(); // different material groups colored with distinct colors (default)

  static void UnitCube();
  static void UnitCubeWireframe();

protected:
  double maxE;
  double maxnu;
  double maxDensity;
  double minE;
  double minnu;
  double minDensity;

  int renderingMode;

  void CubeDeformable(double u0x,double u0y,double u0z,
                      double u1x,double u1y,double u1z,
                      double u2x,double u2y,double u2z,
                      double u3x,double u3y,double u3z,
                      double u4x,double u4y,double u4z,
                      double u5x,double u5y,double u5z,
                      double u6x,double u6y,double u6z,
                      double u7x,double u7y,double u7z
                      );

  void CubeWireframeDeformable(double u0x,double u0y,double u0z,
                        double u1x,double u1y,double u1z,
                        double u2x,double u2y,double u2z,
                        double u3x,double u3y,double u3z,
                        double u4x,double u4y,double u4z,
                        double u5x,double u5y,double u5z,
                        double u6x,double u6y,double u6z,
                        double u7x,double u7y,double u7z
                                        );

   void TetDeformable(double u0x,double u0y,double u0z,
                                         double u1x,double u1y,double u1z,
                                         double u2x,double u2y,double u2z,
                                         double u3x,double u3y,double u3z
                                        );

   void TetWireframeDeformable(double u0x,double u0y,double u0z,
                                         double u1x,double u1y,double u1z,
                                         double u2x,double u2y,double u2z,
                                         double u3x,double u3y,double u3z
                                        );


  void RenderCube(VolumetricMesh * volumetricMesh, int el, int wireframe=0);
  void RenderTet(VolumetricMesh * volumetricMesh, int el, int wireframe=0);

  void JetColorMap(double x, double color[3]);
  void DetermineMaxMin(VolumetricMesh * volumetricMesh);
};

#endif
