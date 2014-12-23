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

#include <float.h>

#ifdef WIN32
  #include <windows.h>
#endif

#include "openGL-headers.h"
#include "renderVolumetricMesh.h"
#include "volumetricMeshENuMaterial.h"
#include "cubicMesh.h"
#include "tetMesh.h"
using namespace std;

// controls how material groups are rendered
#define RENDERVOLUMETRICMESH_RENDERINGMODE_FLAT 0
#define RENDERVOLUMETRICMESH_RENDERINGMODE_DISCRETECOLORS 1
#define RENDERVOLUMETRICMESH_RENDERINGMODE_GRADEDCOLORS 2

RenderVolumetricMesh::RenderVolumetricMesh()
{
  renderingMode = RENDERVOLUMETRICMESH_RENDERINGMODE_DISCRETECOLORS;
}

void RenderVolumetricMesh::DetermineMaxMin(VolumetricMesh * volumetricMesh)
{
  maxE = 0;
  maxnu = 0;
  maxDensity = 0;
  minE = DBL_MAX;
  minnu = DBL_MAX;
  minDensity = DBL_MAX;
  for(int i=0; i < volumetricMesh->getNumElements(); i++)
  {
    // set color based on Young's modulus, Poisson ratio, density
    VolumetricMesh::Material * material = volumetricMesh->getElementMaterial(i);
    double density = material->getDensity();

    VolumetricMesh::ENuMaterial * eNuMaterial =  downcastENuMaterial(material);
    double E;
    double nu;
    if (eNuMaterial == NULL)
    {
      E = 1E6;
      nu = 0.45;
    }
    else
    {
      E = eNuMaterial->getE();
      nu = eNuMaterial->getNu();
    }

    if (E > maxE)
      maxE = E;

    if (nu > maxnu)
      maxnu = nu;

    if(density>maxDensity)
      maxDensity = density;

    if (E < minE)
      minE = E;

    if (nu < minnu)
      minnu = nu;

    if (density<minDensity)
      minDensity = density;
  }

  printf("MaxE: %G MinE: %G\n",maxE,minE);
  printf("Maxnu: %G Minnu: %G\n",maxnu,minnu);
  printf("MaxDensity: %G MinDensity: %G\n",maxDensity,minDensity);
}

void RenderVolumetricMesh::RenderTet(VolumetricMesh * volumetricMesh, int el, int wireframe)
{
  Vec3d v0 = *(volumetricMesh->getVertex(el,0));
  Vec3d v1 = *(volumetricMesh->getVertex(el,1));
  Vec3d v2 = *(volumetricMesh->getVertex(el,2));
  Vec3d v3 = *(volumetricMesh->getVertex(el,3));

  #define RENDERVTX(i) glVertex3f(v##i[0], v##i[1], v##i[2]);
  if (wireframe)
  {
    glBegin(GL_LINES);
      RENDERVTX(0);
      RENDERVTX(1);

      RENDERVTX(1);
      RENDERVTX(2);

      RENDERVTX(2);
      RENDERVTX(0);

      RENDERVTX(0);
      RENDERVTX(3);

      RENDERVTX(1);
      RENDERVTX(3);

      RENDERVTX(2);
      RENDERVTX(3);
    glEnd();
  }
  else
  {
    glBegin(GL_TRIANGLES);
      RENDERVTX(0);
      RENDERVTX(1);
      RENDERVTX(2);

      RENDERVTX(1);
      RENDERVTX(2);
      RENDERVTX(3);

      RENDERVTX(0);
      RENDERVTX(2);
      RENDERVTX(3);

      RENDERVTX(0);
      RENDERVTX(1);
      RENDERVTX(3);
    glEnd();    
  }
}

void RenderVolumetricMesh::RenderCube(VolumetricMesh * volumetricMesh, int el, int wireframe)
{
  // move to vertex 0
  glPushMatrix();
  Vec3d v0 = *(volumetricMesh->getVertex(el,0));
  glTranslated(v0[0],v0[1],v0[2]);

  Vec3d v1 = *(volumetricMesh->getVertex(el,1));
  Vec3d v3 = *(volumetricMesh->getVertex(el,3));
  Vec3d v4 = *(volumetricMesh->getVertex(el,4));
    
  Vec3d axisX = norm(v1-v0);
  Vec3d axisY = norm(v3-v0);
  Vec3d axisZ = norm(v4-v0);

  double M[16] = {axisX[0], axisX[1], axisX[2], 0,
                  axisY[0], axisY[1], axisY[2], 0,
                  axisZ[0], axisZ[1], axisZ[2], 0,
	           0,    0,    0,    1 };

  glMultMatrixd(M);

  double cubeSize = ((CubicMesh*)volumetricMesh)->getCubeSize();
  glScaled(cubeSize, cubeSize, cubeSize);

  if (wireframe)
    UnitCubeWireframe();
  else
  {
    UnitCube();
  }

  glPopMatrix();
}

void RenderVolumetricMesh::JetColorMap(double x, double color[3])
{
  double a; // alpha

  if (x < 0)
  {
    color[0] = 0;
    color[1] = 0;
    color[2] = 0;
    return;
  } 
  else if (x < 0.125) 
  {
    a = x / 0.125;
    color[0] = 0;
    color[1] = 0;
    color[2] = 0.5 + 0.5 * a;
    return;
  }
  else if (x < 0.375) 
  {
    a = (x - 0.125) / 0.25;
    color[0] = 0;
    color[1] = a;
    color[2] = 1;
    return;
  }
  else if (x < 0.625) 
  {         
    a = (x - 0.375) / 0.25;
    color[0] = a;
    color[1] = 1;
    color[2] = 1 - a;
    return;
  }     
  else if (x < 0.875) 
  {
    a = (x - 0.625) / 0.25;
    color[0] = 1;
    color[1] = 1 - a;
    color[2] = 0;
    return;
  }     
  else if (x <= 1.0) 
  {
    a = (x - 0.875) / 0.125;
    color[0] = 1 - 0.5 * a;
    color[1] = 0;
    color[2] = 0;
    return;
  }
  else
  {
    color[0] = 1;
    color[1] = 1;
    color[2] = 1;
    return;
  }
}

void RenderVolumetricMesh::Render(VolumetricMesh * volumetricMesh, int wireframe, double * u)
{
  glDisable(GL_LIGHTING);

  int meshType = 0;
  if (volumetricMesh->getElementType() == CubicMesh::elementType())
    meshType = 1;
  if (volumetricMesh->getElementType() == TetMesh::elementType())
    meshType = 2;

  if (renderingMode == RENDERVOLUMETRICMESH_RENDERINGMODE_DISCRETECOLORS)
  {
    if (volumetricMesh->getNumMaterials() == 0)
    {
      printf("Error: discrete color rendering mode in renderVolumetricMesh called with zero materials.\n");
    }

    map<int,int> actualMaterials;
    for(int ss=0; ss < volumetricMesh->getNumRegions(); ss++)
    {
      int materialIndex = volumetricMesh->getRegion(ss)->getMaterialIndex();
      int setIndex = volumetricMesh->getRegion(ss)->getSetIndex();
      VolumetricMesh::Set * elementSet = volumetricMesh->getSet(setIndex);
      int numElements = elementSet->getNumElements();

      map<int,int> :: iterator iter = actualMaterials.find(materialIndex);
      if (iter == actualMaterials.end())
      {
        // new material
        actualMaterials.insert(make_pair(materialIndex, numElements));
      }
      else
      {
        // existing material
        iter->second += numElements;
      }
    }
    int numActualMaterials = (int)actualMaterials.size();

    multimap<int, int> materialOrderReverse;
    for(map<int, int> :: iterator iter = actualMaterials.begin(); iter != actualMaterials.end(); iter++)
      materialOrderReverse.insert(make_pair(iter->second, iter->first));

    // sort by the number of elements
    map<int,int> materialOrder;
    int counter = 0;
    for(multimap<int, int> :: iterator iter = materialOrderReverse.begin(); iter != materialOrderReverse.end(); iter++)
    {
      materialOrder.insert(make_pair(iter->second, counter));
      counter++;
    }

    double multiplicator = (numActualMaterials == 1 ? 1.0 : 1.0 / (numActualMaterials - 1));
    for(int ss=0; ss < volumetricMesh->getNumRegions(); ss++)
    {
      VolumetricMesh::Region * region = volumetricMesh->getRegion(ss);

      int materialIndex = region->getMaterialIndex(); 
      double color[3];
      //JetColorMap(materialIndex * multiplicator, color); 
      double gray = 1.0 - 0.5 * (numActualMaterials - 1 - materialOrder[materialIndex]) * multiplicator;
      color[0] = gray;
      color[1] = gray;
      color[2] = gray;

      int setIndex = region->getSetIndex(); 
      VolumetricMesh::Set * elementSet = volumetricMesh->getSet(setIndex);
      set<int> elements;
      elementSet->getElements(elements);

      for(set<int> :: iterator iter = elements.begin(); iter != elements.end(); iter++)
      {
        if (wireframe)
          glColor4d(0, 0, 0, 0.8);
        else
          glColor3f(color[0], color[1], color[2]);

        int el = *iter;
   
        if (u == NULL)
        {
          if (meshType == 1)
            RenderCube(volumetricMesh, el, wireframe);

          if (meshType == 2)
            RenderTet(volumetricMesh, el, wireframe);
        }
        else
        {
          #define VER(j) (*volumetricMesh->getVertex(el,j))[0] + u[3*volumetricMesh->getVertexIndex(el, j)+0],\
		         (*volumetricMesh->getVertex(el,j))[1] + u[3*volumetricMesh->getVertexIndex(el, j)+1],\
		         (*volumetricMesh->getVertex(el,j))[2] + u[3*volumetricMesh->getVertexIndex(el, j)+2]
          if (wireframe)
          {
            if (meshType == 1)
              CubeWireframeDeformable(VER(0),VER(1),VER(2),VER(3),
                                      VER(4),VER(5),VER(6),VER(7));

            if (meshType == 2)
              TetWireframeDeformable(VER(0),VER(1),VER(2),VER(3));
          }
          else
          {
            if (meshType == 1)
              CubeDeformable(VER(0),VER(1),VER(2),VER(3),
                             VER(4),VER(5),VER(6),VER(7));

            if (meshType == 2)
              TetDeformable(VER(0),VER(1),VER(2),VER(3));
          }
        }
      }
    }
  }
  else
  {
    for (int i=0; i < volumetricMesh->getNumElements(); i++)
    {
      // set color based on Young's modulus, Poisson ratio, density
      VolumetricMesh::Material * material = volumetricMesh->getElementMaterial(i);
      double density = material->getDensity();

      VolumetricMesh::ENuMaterial * eNuMaterial =  downcastENuMaterial(material);
      double E;
      double nu;
      if (eNuMaterial == NULL)
      {
        E = 1E6;
        nu = 0.45;
      }
      else
      {
        E = eNuMaterial->getE();
        nu = eNuMaterial->getNu();
      }

      double colorR=0.0, colorG=0.0, colorB=0.0;

      if (renderingMode == RENDERVOLUMETRICMESH_RENDERINGMODE_FLAT)
      {
        colorR = colorG = colorB = 1.0;
      }
      else if (renderingMode == RENDERVOLUMETRICMESH_RENDERINGMODE_GRADEDCOLORS)
      {
        if (maxE > minE + 1E-10)
          colorR = (E - minE) / (maxE - minE);
        else
          colorR = 1;
    
        if (maxnu > minnu + 1E-10)
          colorG = (nu - minnu) / (maxnu - minnu);
        else
          colorG = 1;
  
        if (maxDensity > minDensity + 1E-10)
          colorB = (density - minDensity) / (maxDensity - minDensity);
        else
          colorB = 1;
      }
      else
      {
        printf("Error: invalid rendering mode in renderVolumetricMesh!\n");
      }

      if (wireframe)
        glColor4d(0, 0, 0, 0.8);
      else
        glColor3f(colorR, colorG, colorB);

      if (u == NULL)
      {
        if (meshType == 1)
          RenderCube(volumetricMesh, i, wireframe);

        if (meshType == 2)
          RenderTet(volumetricMesh, i, wireframe);
      }
      else
      {
        #define VERA(j) (*volumetricMesh->getVertex(i,j))[0] + u[3*volumetricMesh->getVertexIndex(i, j)+0],\
                        (*volumetricMesh->getVertex(i,j))[1] + u[3*volumetricMesh->getVertexIndex(i, j)+1],\
	                (*volumetricMesh->getVertex(i,j))[2] + u[3*volumetricMesh->getVertexIndex(i, j)+2]

        if (wireframe)
        {
          if (meshType == 1)
            CubeWireframeDeformable(VERA(0),VERA(1),VERA(2),VERA(3),
                           VERA(4),VERA(5),VERA(6),VERA(7));

          if (meshType == 2)
            TetWireframeDeformable(VERA(0),VERA(1),VERA(2),VERA(3));
        }
        else
        {
          if (meshType == 1)
            CubeDeformable(VERA(0),VERA(1),VERA(2),VERA(3),
                           VERA(4),VERA(5),VERA(6),VERA(7));

          if (meshType == 2)
            TetDeformable(VERA(0),VERA(1),VERA(2),VERA(3));
        }
      }
    }
  }
}

void RenderVolumetricMesh::RenderWireframe(VolumetricMesh * volumetricMesh)
{
  Render(volumetricMesh, 1, NULL);
}

void RenderVolumetricMesh::RenderDeformation(VolumetricMesh * volumetricMesh, double * u)
{
  Render(volumetricMesh, 0, u);
}

void RenderVolumetricMesh::RenderWireframeDeformation(VolumetricMesh * volumetricMesh, double * u)
{
  Render(volumetricMesh, 1, u);
}

void RenderVolumetricMesh::RenderVertexDeformed(VolumetricMesh * volumetricMesh, int ver, double * U)
{
  glBegin(GL_POINTS);
    glVertex3f((*volumetricMesh->getVertex(ver))[0]+U[3*(ver)+0],
	       (*volumetricMesh->getVertex(ver))[1]+U[3*(ver)+1],
               (*volumetricMesh->getVertex(ver))[2]+U[3*(ver)+2]);
  glEnd();
}


void RenderVolumetricMesh::RenderVertices(VolumetricMesh * volumetricMesh)
{
  int i;
  glBegin(GL_POINTS);
  for (i=0; i < volumetricMesh->getNumVertices(); i++)
  {
     Vec3d * vertex = volumetricMesh->getVertex(i);
     glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
  }
  glEnd();
}

void RenderVolumetricMesh::RenderVertices
  (VolumetricMesh * volumetricMesh, std::set<int> * vertices, bool oneIndexed)
{
  glBegin(GL_POINTS);
  int offset = (oneIndexed ? -1 : 0);
  std::set<int> :: iterator iter;
  for(iter = vertices->begin(); iter != vertices->end(); iter++)
  {
    Vec3d * vertex = volumetricMesh->getVertex(*iter + offset);
    glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
  }
  glEnd();
}

void RenderVolumetricMesh::RenderVertices(VolumetricMesh * volumetricMesh, 
                                  int * vertices, int numVertices, bool oneIndexed)
{
  glBegin(GL_POINTS);
  int offset = (oneIndexed ? -1 : 0);
  for (int i=0; i < numVertices; i++)
  {
    Vec3d * vertex = volumetricMesh->getVertex(vertices[i] + offset);
    glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
  }
  glEnd();
}

void RenderVolumetricMesh::SelectRenderVertices(VolumetricMesh * volumetricMesh)
{
  for (int i=0; i < volumetricMesh->getNumVertices(); i++)
  {
    glLoadName(i+1);
    glBegin(GL_POINTS);
      Vec3d * vertex = volumetricMesh->getVertex(i);
      glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
    glEnd();
  }
}

void RenderVolumetricMesh::DrawSelectedPoints(VolumetricMesh * volumetricMesh, int * selectedVertices, int numSelectedVertices)
{
  glBegin(GL_POINTS);
  for (int i=0; i < numSelectedVertices; i++)
  {
    Vec3d * vertex = volumetricMesh->getVertex(selectedVertices[i]-1);
    glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
  }
  glEnd();
}

void RenderVolumetricMesh::DrawUnselectedPoints(VolumetricMesh * volumetricMesh, int * selectionArray)
{
  glBegin(GL_POINTS);
  for (int i=0; i < volumetricMesh->getNumVertices(); i++)
  {
    if (selectionArray[i+1] != 0)
      continue;
    Vec3d * vertex = volumetricMesh->getVertex(i);
    glVertex3f((*vertex)[0],(*vertex)[1],(*vertex)[2]);
  }
  glEnd();
}

void RenderVolumetricMesh::RenderSolidAndWireframe(VolumetricMesh * volumetricMesh)
{
  // do the polygon offset trick

  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (2., 1.);
  Render(volumetricMesh);
  glDisable (GL_POLYGON_OFFSET_FILL);

  glLineWidth(1);
  RenderWireframe(volumetricMesh);
}

void RenderVolumetricMesh::SetFlatRenderingMode()
{
  renderingMode = RENDERVOLUMETRICMESH_RENDERINGMODE_FLAT;
}

void RenderVolumetricMesh::SetGradedRenderingMode(VolumetricMesh * volumetricMesh)
{
  renderingMode = RENDERVOLUMETRICMESH_RENDERINGMODE_GRADEDCOLORS;
  DetermineMaxMin(volumetricMesh);
}

void RenderVolumetricMesh::SetDiscreteRenderingMode()
{
  renderingMode = RENDERVOLUMETRICMESH_RENDERINGMODE_DISCRETECOLORS;
}

void RenderVolumetricMesh::RenderSolidAndWireframeDeformation(VolumetricMesh * volumetricMesh, double * U)
{
  // do the polygon offset trick
  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (2., 1.);
  RenderDeformation(volumetricMesh,U);
  glDisable (GL_POLYGON_OFFSET_FILL);

  glLineWidth(1);
  RenderWireframeDeformation(volumetricMesh,U);
}

void RenderVolumetricMesh::RenderVertexLabels(VolumetricMesh * volumetricMesh)
{
  RenderVertexLabels(volumetricMesh,0, volumetricMesh->getNumVertices());
}

void RenderVolumetricMesh::RenderVertexLabels(VolumetricMesh * volumetricMesh, int start, int end)
{
  // show point labels
  // labels are printed out in the range 1... , not 0...
  for (int i=start; i< end; i++)
  {
    Vec3d * vertex = volumetricMesh->getVertex(i);
    print_bitmap_integer((*vertex)[0],(*vertex)[1],(*vertex)[2],i+1);
  }
}

void RenderVolumetricMesh::UnitCube()
{
  glBegin(GL_TRIANGLES);

  glNormal3f(0,-1,0);  

  glVertex3f(0,0,0); // front
  glVertex3f(1,0,0);
  glVertex3f(1,0,1);

  glVertex3f(0,0,0);
  glVertex3f(1,0,1);
  glVertex3f(0,0,1);

  glNormal3f(0,1,0);  

  glVertex3f(0,1,0); // back
  glVertex3f(1,1,1);
  glVertex3f(1,1,0);

  glVertex3f(0,1,1);
  glVertex3f(1,1,1);
  glVertex3f(0,1,0);

  glNormal3f(1,0,0);  

  glVertex3f(1,0,0); // right
  glVertex3f(1,1,0);
  glVertex3f(1,1,1);

  glVertex3f(1,0,0);
  glVertex3f(1,1,1);
  glVertex3f(1,0,1);

  glNormal3f(-1,0,0);  

  glVertex3f(0,0,0); // left
  glVertex3f(0,1,1);
  glVertex3f(0,1,0);

  glVertex3f(0,0,1);
  glVertex3f(0,1,1);
  glVertex3f(0,0,0);

  glNormal3f(0,0,1);  

  glVertex3f(0,0,1); // top
  glVertex3f(1,0,1);
  glVertex3f(1,1,1);

  glVertex3f(0,0,1);
  glVertex3f(1,1,1);
  glVertex3f(0,1,1);

  glNormal3f(0,0,-1);  

  glVertex3f(0,0,0); // bottom
  glVertex3f(1,1,0);
  glVertex3f(1,0,0);

  glVertex3f(0,1,0);
  glVertex3f(1,1,0);
  glVertex3f(0,0,0);

  glEnd();
}

void RenderVolumetricMesh::UnitCubeWireframe()
{
  glBegin(GL_LINES); 
    glVertex3f(0,0,0);
    glVertex3f(1,0,0);
    glVertex3f(0,1,0);
    glVertex3f(1,1,0);
    glVertex3f(0,0,0);
    glVertex3f(0,1,0);
    glVertex3f(1,0,0);
    glVertex3f(1,1,0);

    glVertex3f(0,0,1);
    glVertex3f(1,0,1);
    glVertex3f(0,1,1);
    glVertex3f(1,1,1);
    glVertex3f(0,0,1);
    glVertex3f(0,1,1);
    glVertex3f(1,0,1);
    glVertex3f(1,1,1);

    glVertex3f(0,0,0);
    glVertex3f(0,0,1);

    glVertex3f(0,1,0);
    glVertex3f(0,1,1);

    glVertex3f(1,0,0);
    glVertex3f(1,0,1);

    glVertex3f(1,1,0);
    glVertex3f(1,1,1);
  glEnd();
}

void RenderVolumetricMesh::CubeDeformable(double u0x,double u0y,double u0z,
					double u1x,double u1y,double u1z,
					double u2x,double u2y,double u2z,
					double u3x,double u3y,double u3z,
					double u4x,double u4y,double u4z,
					double u5x,double u5y,double u5z,
					double u6x,double u6y,double u6z,
					double u7x,double u7y,double u7z
					)
{
  glBegin(GL_TRIANGLES);

  glNormal3f(0,-1,0);  

  glVertex3f(u0x,u0y,u0z); // front
  glVertex3f(u1x,u1y,u1z);
  glVertex3f(u5x,u5y,u5z);

  glVertex3f(u0x,u0y,u0z);
  glVertex3f(u5x,u5y,u5z);
  glVertex3f(u4x,u4y,u4z);

  glNormal3f(u3x,u3y,u3z);  

  glVertex3f(u3x,u3y,u3z); // back
  glVertex3f(u6x,u6y,u6z);
  glVertex3f(u2x,u2y,u2z);

  glVertex3f(u7x,u7y,u7z);
  glVertex3f(u6x,u6y,u6z);
  glVertex3f(u3x,u3y,u3z);

  glNormal3f(u1x,u1y,u1z);  

  glVertex3f(u1x,u1y,u1z); // right
  glVertex3f(u2x,u2y,u2z);
  glVertex3f(u6x,u6y,u6z);

  glVertex3f(u1x,u1y,u1z);
  glVertex3f(u6x,u6y,u6z);
  glVertex3f(u5x,u5y,u5z);

  glNormal3f(-1,0,0);  

  glVertex3f(u0x,u0y,u0z); // left
  glVertex3f(u7x,u7y,u7z);
  glVertex3f(u3x,u3y,u3z);

  glVertex3f(u4x,u4y,u4z);
  glVertex3f(u7x,u7y,u7z);
  glVertex3f(u0x,u0y,u0z);

  glNormal3f(u4x,u4y,u4z);  

  glVertex3f(u4x,u4y,u4z); // top
  glVertex3f(u5x,u5y,u5z);
  glVertex3f(u6x,u6y,u6z);

  glVertex3f(u4x,u4y,u4z);
  glVertex3f(u6x,u6y,u6z);
  glVertex3f(u7x,u7y,u7z);

  glNormal3f(0,0,-1);  

  glVertex3f(u0x,u0y,u0z); // bottom
  glVertex3f(u2x,u2y,u2z);
  glVertex3f(u1x,u1y,u1z);

  glVertex3f(u3x,u3y,u3z);
  glVertex3f(u2x,u2y,u2z);
  glVertex3f(u0x,u0y,u0z);

  glEnd();
}

void RenderVolumetricMesh::CubeWireframeDeformable(double u0x,double u0y,double u0z,
    			double u1x,double u1y,double u1z,
			double u2x,double u2y,double u2z,
			double u3x,double u3y,double u3z,
			double u4x,double u4y,double u4z,
			double u5x,double u5y,double u5z,
			double u6x,double u6y,double u6z,
			double u7x,double u7y,double u7z
					)
{
  glBegin(GL_LINES); 
    glVertex3f(u0x,u0y,u0z);
    glVertex3f(u1x,u1y,u1z);
    glVertex3f(u3x,u3y,u3z);
    glVertex3f(u2x,u2y,u2z);
    glVertex3f(u0x,u0y,u0z);
    glVertex3f(u3x,u3y,u3z);
    glVertex3f(u1x,u1y,u1z);
    glVertex3f(u2x,u2y,u2z);

    glVertex3f(u4x,u4y,u4z);
    glVertex3f(u5x,u5y,u5z);
    glVertex3f(u7x,u7y,u7z);
    glVertex3f(u6x,u6y,u6z);
    glVertex3f(u4x,u4y,u4z);
    glVertex3f(u7x,u7y,u7z);
    glVertex3f(u5x,u5y,u5z);
    glVertex3f(u6x,u6y,u6z);

    glVertex3f(u0x,u0y,u0z);
    glVertex3f(u4x,u4y,u4z);

    glVertex3f(u3x,u3y,u3z);
    glVertex3f(u7x,u7y,u7z);

    glVertex3f(u1x,u1y,u1z);
    glVertex3f(u5x,u5y,u5z);

    glVertex3f(u2x,u2y,u2z);
    glVertex3f(u6x,u6y,u6z);
  glEnd();
}

void RenderVolumetricMesh::TetDeformable(double u0x,double u0y,double u0z,
					 double u1x,double u1y,double u1z,
					 double u2x,double u2y,double u2z,
					 double u3x,double u3y,double u3z
					)
{
  #define RENDERVTXALT(i) glVertex3f(u##i##x, u##i##y, u##i##z);
  glBegin(GL_TRIANGLES);
    RENDERVTXALT(0);
    RENDERVTXALT(1);
    RENDERVTXALT(2);

    RENDERVTXALT(1);
    RENDERVTXALT(2);
    RENDERVTXALT(3);

    RENDERVTXALT(0);
    RENDERVTXALT(2);
    RENDERVTXALT(3);

    RENDERVTXALT(0);
    RENDERVTXALT(1);
    RENDERVTXALT(3);
  glEnd();    
}

void RenderVolumetricMesh::TetWireframeDeformable(double u0x,double u0y,double u0z,
					 double u1x,double u1y,double u1z,
					 double u2x,double u2y,double u2z,
					 double u3x,double u3y,double u3z
					)
{
  glBegin(GL_LINES);
    RENDERVTXALT(0);
    RENDERVTXALT(1);

    RENDERVTXALT(1);
    RENDERVTXALT(2);

    RENDERVTXALT(2);
    RENDERVTXALT(0);

    RENDERVTXALT(0);
    RENDERVTXALT(3);

    RENDERVTXALT(1);
    RENDERVTXALT(3);

    RENDERVTXALT(2);
    RENDERVTXALT(3);
  glEnd();
}

