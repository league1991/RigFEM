/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2014 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
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

#ifdef WIN32
  #include <windows.h>
#endif
#include "openGL-headers.h"
#include "renderSprings.h"

void RenderSprings::Render(MassSpringSystem * massSpringSystem, double * u)
{
  for(int edge=0; edge < massSpringSystem->numEdges; edge++)
  {
    int group = massSpringSystem->edgeGroups[edge]; 
    
    // get color
    double color[3];
    double lookupColorIndex = 1.0 * group / (massSpringSystem->numMaterialGroups + 1) ;
    JetColorMap(lookupColorIndex, color);

    int vtxA = massSpringSystem->edges[2*edge+0];
    int vtxB = massSpringSystem->edges[2*edge+1];
    
    double posA[3] = { massSpringSystem->restPositions[3*vtxA+0] + u[3*vtxA+0], massSpringSystem->restPositions[3*vtxA+1] + u[3*vtxA+1], massSpringSystem->restPositions[3*vtxA+2] + u[3*vtxA+2] };
    double posB[3] = { massSpringSystem->restPositions[3*vtxB+0] + u[3*vtxB+0], massSpringSystem->restPositions[3*vtxB+1] + u[3*vtxB+1], massSpringSystem->restPositions[3*vtxB+2] + u[3*vtxB+2] };

    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
      glVertex3f(posA[0], posA[1], posA[2]);
      glVertex3f(posB[0], posB[1], posB[2]);
    glEnd();
  }
}

void RenderSprings::JetColorMap(double x, double color[3])
{
  double a;

  if(x < 0)
  {
    color[0] = color[1] = color[2] = 0.0;
    return;
  }
  else if (x < 0.125) 
  {
    a = x/0.125;
    color[0] = color[1] = 0;
    color[2] = 0.5 + 0.5 * a;
    return;
  }
  else if (x < 0.375) 
  {
    a = (x - 0.125)/0.25;
    color[0] = 0;
    color[1] = a;
    color[2] = 1;
    return;
  }
  else if (x < 0.625) 
  {
    a = (x - 0.375)/0.25;
    color[0] = a;
    color[1] = 1;
    color[2] = 1 - a;
    return;
  }
  else if (x < 0.875) 
  {
    a = (x - 0.625)/0.25;
    color[0] = 1;
    color[1] = 1 - a;
    color[2] = 0;
    return;
  }
  else if (x <= 1.0) 
  {
    a = (x - 0.875)/0.125;
    color[0] = 1 - 0.5 * a;
    color[1] = color[2] = 0.0;
    return;
  }
  else 
  {
    color[0] = color[1] = color[2] = 1.0;
    return;
  }
}

