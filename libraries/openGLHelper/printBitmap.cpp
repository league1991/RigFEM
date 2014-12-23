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

#include "printBitmap.h"
#include "openGLHelper.h"
#ifdef WIN32
//  #include "windows.h"
#endif

#include "openGL-headers.h"

void print_bitmap_string(float x, float y, float z, const char * s)
{
  glRasterPos3f(x,y,z);
  if (s && strlen(s)) 
  {
    while (*s) 
    {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *s);
      s++;
    }
  }
}

void print_bitmap_integer(float x,float y, float z, int i)
{
  char s[50];
  sprintf(s,"%d",i);
  print_bitmap_string(x,y,z,s);
}

