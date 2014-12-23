/*
  * Copyright (c) 2007, Carnegie Mellon University
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of Carnegie Mellon University, nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "configFile.h"
#include "lighting.h"
#ifdef WIN32
  #include <windows.h>
#endif
#include "openGL-headers.h"

Lighting::Lighting(const char * configFilename)
{
  if (glGetError() != GL_NO_ERROR)
    printf("Warning: error in the OpenGL state at the beginning of lighting constructor.\n");

  ConfigFile configFile;

  configFile.addOptionOptional("globalAmbientIntensity", &ambientIntensity, 0.0f);
  configFile.addOptionOptional("localViewer", &localViewer, true);
  configFile.addOptionOptional("twoSidedLighting", &twoSidedLighting, false);
  configFile.addOptionOptional("enableAmbientTerm", &enableAmbientTerm, true);
  configFile.addOptionOptional("enableDiffuseTerm", &enableDiffuseTerm, true);
  configFile.addOptionOptional("enableSpecularTerm", &enableSpecularTerm, true);

  for(int light=0; light<8; light++)
  {
    char lightCh = '0' + light;
    char optionName[128];

    sprintf(optionName, "lightEnabled_%c", lightCh);
    configFile.addOptionOptional(optionName, &lightEnabled[light], false);

    sprintf(optionName, "position_%c_X", lightCh);
    configFile.addOptionOptional(optionName, &lightPosition[4*light+0], 0.0f);
    sprintf(optionName, "position_%c_Y", lightCh);
    configFile.addOptionOptional(optionName, &lightPosition[4*light+1], 0.0f);
    sprintf(optionName, "position_%c_Z", lightCh);
    configFile.addOptionOptional(optionName, &lightPosition[4*light+2], 0.0f);
    sprintf(optionName, "position_%c_DIR", lightCh);
    configFile.addOptionOptional(optionName, &lightPosition[4*light+3], 1.0f);

    sprintf(optionName, "lightIntensity_%c", lightCh);
    configFile.addOptionOptional(optionName, &lightIntensity[light], 1.0f);
  }
  
  if (configFile.parseOptions(configFilename) != 0)
    throw 1;

  //configFile.printOptions();

  #define TURNONLIGHT(i)\
   if(lightEnabled[i])\
     glEnable(GL_LIGHT##i);\
   else\
     glDisable(GL_LIGHT##i);\

  TURNONLIGHT(0);
  TURNONLIGHT(1);
  TURNONLIGHT(2);
  TURNONLIGHT(3);
  TURNONLIGHT(4);
  TURNONLIGHT(5);
  TURNONLIGHT(6);
  TURNONLIGHT(7);

  if (glGetError() != GL_NO_ERROR)
    printf("Warning: error in the OpenGL state in the Lighting constructor.\n");

  for(int light=0; light<8; light++)
  {
    if(!lightEnabled[light])
      continue;

    if (enableAmbientTerm) 
    {
      lKa[light*4+0] = lKa[light*4+1] = lKa[light*4+2] = lightIntensity[light];
    }
    else
    {
      lKa[light*4+0] = lKa[light*4+1] = lKa[light*4+2] = 0.0f;
    }

    if (enableDiffuseTerm)
    {
      lKd[light*4+0] = lKd[light*4+1] = lKd[light*4+2] = lightIntensity[light];
    }
    else
    {
      lKd[light*4+0] = lKd[light*4+1] = lKd[light*4+2] = 0.0f;
    }

    if (enableSpecularTerm)
    {
      lKs[light*4+0] = lKs[light*4+1] = lKs[light*4+2] = lightIntensity[light];
    }
    else
    {
      lKs[light*4+0] = lKs[light*4+1] = lKs[light*4+2] = 0.0f;
    }

    // set alpha to 1.0
    lKa[light*4+3] = lKd[light*4+3] = lKs[light*4+3] = 1.0;
  }

  GLfloat aGa[] = { ambientIntensity, ambientIntensity, ambientIntensity, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, localViewer);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, twoSidedLighting);
}

void Lighting::LightScene()
{
  #define LIGHTSETUP(i)\
  if (lightEnabled[i])\
  {\
    glLightfv(GL_LIGHT##i, GL_POSITION, &lightPosition[4*i]);\
    glLightfv(GL_LIGHT##i, GL_AMBIENT, &lKa[4*i]);\
    glLightfv(GL_LIGHT##i, GL_DIFFUSE, &lKd[4*i]);\
    glLightfv(GL_LIGHT##i, GL_SPECULAR, &lKs[4*i]);\
  }

  LIGHTSETUP(0);
  LIGHTSETUP(1);
  LIGHTSETUP(2);
  LIGHTSETUP(3);
  LIGHTSETUP(4);
  LIGHTSETUP(5);
  LIGHTSETUP(6);
  LIGHTSETUP(7);

  glEnable(GL_LIGHTING);
}

