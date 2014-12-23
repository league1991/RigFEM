/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "openGLHelper" library , Copyright (C) 2007 CMU, 2009 MIT, 2014 USC   *
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

#include "openGLHelper.h"

void *font = GLUT_BITMAP_9_BY_15;
void *fonts[] =
{
  GLUT_BITMAP_9_BY_15,
  GLUT_BITMAP_TIMES_ROMAN_10,
  GLUT_BITMAP_TIMES_ROMAN_24
};


void UnitCube()
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

void UnitCubeWireframe()
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

void RenderWireframeBox(double bmin[3], double bmax[3])
{
  glBegin(GL_LINES); 
    glVertex3f(bmin[0],bmin[1],bmin[2]);
    glVertex3f(bmax[0],bmin[1],bmin[2]);
    glVertex3f(bmin[0],bmax[1],bmin[2]);
    glVertex3f(bmax[0],bmax[1],bmin[2]);
    glVertex3f(bmin[0],bmin[1],bmin[2]);
    glVertex3f(bmin[0],bmax[1],bmin[2]);
    glVertex3f(bmax[0],bmin[1],bmin[2]);
    glVertex3f(bmax[0],bmax[1],bmin[2]);

    glVertex3f(bmin[0],bmin[1],bmax[2]);
    glVertex3f(bmax[0],bmin[1],bmax[2]);
    glVertex3f(bmin[0],bmax[1],bmax[2]);
    glVertex3f(bmax[0],bmax[1],bmax[2]);
    glVertex3f(bmin[0],bmin[1],bmax[2]);
    glVertex3f(bmin[0],bmax[1],bmax[2]);
    glVertex3f(bmax[0],bmin[1],bmax[2]);
    glVertex3f(bmax[0],bmax[1],bmax[2]);

    glVertex3f(bmin[0],bmin[1],bmin[2]);
    glVertex3f(bmin[0],bmin[1],bmax[2]);

    glVertex3f(bmin[0],bmax[1],bmin[2]);
    glVertex3f(bmin[0],bmax[1],bmax[2]);

    glVertex3f(bmax[0],bmin[1],bmin[2]);
    glVertex3f(bmax[0],bmin[1],bmax[2]);

    glVertex3f(bmax[0],bmax[1],bmin[2]);
    glVertex3f(bmax[0],bmax[1],bmax[2]);
  glEnd();

}

void RenderSphere(float x, float y, float z)
{
  glPushMatrix();
  glTranslatef(x,y,z);
  glutSolidSphere(0.02,20,20);
  glPopMatrix();
}

void OutputText(int x, int y, const char *string)
{
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++) 
  {
    glutBitmapCharacter(font, string[i]);
  }
}

void PrintGLerror( const char *msg )
{
 GLenum errCode;
 const GLubyte *errStr;

 if ((errCode = glGetError()) != GL_NO_ERROR) 
 {
    errStr = gluErrorString(errCode);
    printf("OpenGL ERROR: %s: %s\n", errStr, msg);
	exit(1);
 }
}



/*
void RenderTexturesIntoDisplayList()

{
  glPushMatrix();
  //glTranslatef(xLow, yLow, zLow);
  //glScalef(xSize, ySize, zSize);

  // draw the walls
  glEnable(GL_TEXTURE_2D);
  
  glBindTexture(GL_TEXTURE_2D,texName3);

  glBegin(GL_QUADS); // left wall
  glTexCoord2f(0,0); glVertex3f(0 ,0 , 0 );
  glTexCoord2f(1,0); glVertex3f(0 ,1 , 0 );
  glTexCoord2f(1,1); glVertex3f(0 ,1 , 1 );
  glTexCoord2f(0,1); glVertex3f(0 ,0 , 1 );
  glEnd();

  glBegin(GL_QUADS); // far wall
  glTexCoord2f(0,0); glVertex3f(0 ,1 ,0 );
  glTexCoord2f(1,0); glVertex3f(1 ,1 ,0 );
  glTexCoord2f(1,1); glVertex3f(1 ,1 ,1 );
  glTexCoord2f(0,1); glVertex3f(0 ,1 ,1 );
  glEnd();
  
  glBegin(GL_QUADS); // right wall
  glTexCoord2f(0,0); glVertex3f(1 ,1 , 0 );
  glTexCoord2f(1,0); glVertex3f(1 ,0 , 0 );
  glTexCoord2f(1,1); glVertex3f(1 ,0 , 1 );
  glTexCoord2f(0,1); glVertex3f(1 ,1 , 1 );
  glEnd();

  glBegin(GL_QUADS); // near wall
  glTexCoord2f(0,0); glVertex3f(1 ,0 ,0 );
  glTexCoord2f(1,0); glVertex3f(0 ,0 ,0 );
  glTexCoord2f(1,1); glVertex3f(0 ,0 ,1 );
  glTexCoord2f(0,1); glVertex3f(1 ,0 ,1 );
  glEnd();

  glBindTexture(GL_TEXTURE_2D,texName2);
  
  glBegin(GL_QUADS); // floor
  glTexCoord2f(0,0); glVertex3f(0 ,0 ,0 );
  glTexCoord2f(1,0); glVertex3f(1 ,0 ,0 );
  glTexCoord2f(1,1); glVertex3f(1 ,1 ,0 );
  glTexCoord2f(0,1); glVertex3f(0 ,1 ,0 );
  glEnd();

  glBindTexture(GL_TEXTURE_2D,texName1);

  glBegin(GL_QUADS); // ceiling
  glTexCoord2f(0,0); glVertex3f(0 ,1 ,1 );
  glTexCoord2f(1,0); glVertex3f(1 ,1 ,1 );
  glTexCoord2f(1,1); glVertex3f(1 ,0 ,1 );
  glTexCoord2f(0,1); glVertex3f(0 ,0 ,1 );

  glEnd();


  glDisable(GL_TEXTURE_2D);

  glPopMatrix();

}
*/


void DetermineCameraParameters(double centerX, double centerY, double centerZ, double modelRadius,
							   double * focusX, double * focusY, double * focusZ,
							   double * cameraRadius, double * zNear, double * zFar)
{
  *focusX = centerX;
  *focusY = centerY;
  *focusZ = centerZ;

  *cameraRadius = modelRadius * 3;
  
  *zNear = *cameraRadius * 0.01;
  *zFar = *cameraRadius * 100;
}


void DrawArrow( float px, float py, float pz,
    float nx, float ny, float nz,
    double arrowEndWidth, double arrowEndLength )
{
  GLdouble normal[3], cross[3], zaxis[3];
  GLdouble angle;

  #define DegreesToRadians     (3.14159265 / (GLfloat) 180.0)

  normal[0] = nx; normal[1] = ny; normal[2] = nz;
  double len = sqrt(nx*nx + ny*ny + nz*nz);

  if (len < 1E-6)
	return;


  glPushMatrix();
    glTranslatef( px, py, pz );
    glBegin( GL_LINES );
      glVertex3f( 0.0, 0.0, 0.0 );
      glVertex3f( nx, ny, nz );
    glEnd();

    // normalize the normal vector 
	normal[0] /= len;
	normal[1] /= len;
	normal[2] /= len;

    // determine angle between z axis and normal 
    zaxis[0] = 0; zaxis[1] = 0; zaxis[2] = 1;
    angle = acos( zaxis[0]*normal[0] + zaxis[1]*normal[1] + zaxis[2]*normal[2] )/DegreesToRadians;
         
    if ( angle != 0.0 ) 
	{
      // find the axis of rotation 
      CROSSPRODUCT( zaxis[0], zaxis[1], zaxis[2], 
		            normal[0], normal[1], normal[2], 
					cross[0], cross[1], cross[2] );
      glRotatef( angle, cross[0], cross[1], cross[2] );
    }

    // move to end of normal vector 
    glTranslatef( 0.0, 0.0, len );
	/*
	glScalef(0.3,0.3,0.3);
    #ifdef EIFFEL140VOX
  	  glScalef(2.0,2.0,2.0);
    #endif
    #ifdef SPOON100VOX
  	  glScalef(2.0,2.0,2.0);
    #endif
	*/
    glutSolidCone( len*arrowEndWidth, len*arrowEndLength, 12, 1 );
  glPopMatrix();
}

void JetColorMap(double x, double color[3])
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

void TransparentSphere(GLuint solidSphereList, GLuint wireSphereList, double x, double y, double z, double radius)
{
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glPushMatrix();
  glTranslatef(x,y,z);
  glColor3f(1,0,0);
  glBegin(GL_POINTS);
    glVertex3f(0,0,0);
  glEnd();
  //glColor4f(1,1,1,0.35);
  glScalef(radius,radius,radius);
  glCallList(solidSphereList);
  //glutSolidSphere(radius,25,25);
  glDisable(GL_BLEND);
  //glColor3f(0,0,0);
  glCallList(wireSphereList);
  //glutWireSphere(radius,25,25);
  glPopMatrix();
}

void BuildSphereDisplayList(GLuint * solidSphereList, GLuint * wireSphereList)
{ 
  *solidSphereList = glGenLists(1);
  glNewList(*solidSphereList, GL_COMPILE);
    glColor4f(1,1,1,0.35);
    glutSolidSphere(1.0,25,25);
  glEndList();

  *wireSphereList = glGenLists(1);
  glNewList(*wireSphereList, GL_COMPILE);
    glColor3f(0,0,0);
    glutWireSphere(1.0,25,25);
  glEndList();
} 

char * DuplicateString(const char * s)
{
  // strdup sometimes causes problems, so we resort to this
  char * p = (char*) malloc (sizeof(char) * (strlen(s) + 1));
  memcpy(p, s, sizeof(char) * (strlen(s) + 1));
  return p;
}

