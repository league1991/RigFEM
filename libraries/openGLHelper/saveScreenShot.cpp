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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "openGLHelper.h"
#include "openGL-headers.h"
#include "saveScreenShot.h"

/* Write a screenshot to the specified filename */

void Screenshot::SaveScreenshot(const char * filename, ImageIO::fileFormatType fileFormat, int windowWidth, int windowHeight)
{
  if (filename == NULL)
    return;

  printf("Saving screenshot to: %s\n", filename);

  unsigned char * pixels = (unsigned char*) malloc (sizeof(unsigned char) * 3 * windowWidth * windowHeight);

  for (int i=0; i<windowHeight; i++)
    glReadPixels(0, i, windowWidth, 1, GL_RGB, GL_UNSIGNED_BYTE, &pixels[3 * windowWidth * i]);

  ImageIO imageIO(windowWidth, windowHeight, 3, pixels);
  ImageIO::errorType errorCode = imageIO.save(filename, fileFormat);

  if (errorCode == ImageIO::OK)
    printf("File saved successfully.\n");
  else
    printf("Error in saving.\n");

  free(pixels);
}

void Screenshot::SaveStencilBuffer(const char *filename, ImageIO::fileFormatType fileFormat, int windowWidth, int windowHeight, int rescale)
{
  if (filename == NULL)
    return;

  // Allocate a picture buffer
  unsigned char * buffer = (unsigned char*) calloc (windowWidth * windowHeight, sizeof(unsigned char));

  printf("Saving stencil buffer screenshot to: %s\n", filename);

  //for (int i=windowHeight-1; i>=0; i--) 
    //glReadPixels(0, windowHeight-1-i, windowWidth, 1, GL_RGB, GL_UNSIGNED_BYTE, &buffer[i*windowWidth*3] );

  //glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  //glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, windowWidth, windowHeight, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, buffer);

  if (rescale)
  {
    int maxValue = 0;
    for(int i=0; i<windowWidth * windowHeight; i++)
    {
      if (buffer[i] > maxValue)
        maxValue = buffer[i];
    }
    if (maxValue > 0)
    {
      for(int i=0; i<windowWidth * windowHeight; i++)
        buffer[i] = (int)(1.0 * buffer[i] / maxValue * 255.0);
    }
    else
    {
      printf("Warning: max stencil buffer value is 0.\n");
    }
  }

  unsigned char * buffer3 = (unsigned char*) calloc (windowWidth * windowHeight, 3 * sizeof(unsigned char));
  for(int i=0; i<windowWidth * windowHeight; i++)
  {
    buffer3[3*i+0] = buffer[i];
    buffer3[3*i+1] = buffer[i];
    buffer3[3*i+2] = buffer[i];
  }
  free(buffer);

  ImageIO imageIO(windowWidth, windowHeight, 3, buffer3);
  ImageIO::errorType errorCode = imageIO.save(filename, fileFormat);

  if (errorCode == ImageIO::OK)
    printf("File saved successfully.\n");
  else
    printf("Error in saving.\n");

  free(buffer3);
}

