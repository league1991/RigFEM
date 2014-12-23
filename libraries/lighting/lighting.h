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

  Code author: Jernej Barbic
  Research: Jernej Barbic, Doug L. James
  Version: 1.0

  A class that makes it possible to read OpenGL lighting parameters from a text configuration file.
  Usage is simple: read the configuration file during initialization (using the constructor), 
  then call "LightScene" inside your OpenGL display routine, after setting up the modelview and projection matrices.
  
  The configuration file syntax supports the following parameters (listed with their types and default values):

  Global parameters:
  ------------------

  globalAmbientIntensity (float, 0.0)
  localViewer (bool, true)
  twoSidedLighting (bool, false)
  enableAmbientTerm (bool, true)
  enableDiffuseTerm (bool, true)
  enableSpecularTerm (bool, true)
  
  Light-specific parameters:
  --------------------------

  Up to 8 simultaneous lights are supported. 
  Each light has the following parameters (listed here for light 0) :
  lightEnabled_0 (bool, false)
  position_0_X (float, 0.0)
  position_0_Y (float, 0.0)
  position_0_Z (float, 0.0)
  position_0_DIR (float, 1.0)
  lightIntensity_0 (float, 1.0)

  Light 0 is located at (position_0_X, position_0_Y, position_0_Z).
  position_0_DIR should be 0 for directional light, and 1 for a point light source

  For the other lights, replace "0" with one of "1", "2", ..., "7".
  For example: lightEnabled_4 controls whether light 4 is enabled or not.

  See also example.lighting for an example of a lighting configuration file.
*/

#ifndef _LIGHTING_H_
#define _LIGHTING_H_

#include <stdlib.h>
#include <string.h>
#include "macros.h"

class Lighting
{
public:

  // read OpenGL lighting parameters from a configuration file
  Lighting(const char * configFilename);

  // call this inside your OpenGL display routine, after setting up the modelview and projection matrices
  void LightScene(); 

  inline void GetLightPosition(int lightID, float * pos);
  inline void GetLightAmbientTerm(int lightID, float * Ka);
  inline void GetLightDiffuseTerm(int lightID, float * Kd);
  inline void GetLightSpecularTerm(int lightID, float * Ks);
  inline float GetLightIntensity(int lightID);
  inline float GetAmbientIntensity();
  inline bool IsLightEnabled(int lightID);
  inline bool IsAmbientEnabled() { return enableAmbientTerm; }
  inline bool IsDiffuseEnabled() { return enableDiffuseTerm; }
  inline bool IsSpecularEnabled() { return enableSpecularTerm; }

protected:
  float ambientIntensity;
  bool localViewer, twoSidedLighting;
  bool lightEnabled[8];
  float lightPosition[32];
  float lKa[32], lKd[32], lKs[32];
  float lightIntensity[8];
  bool enableAmbientTerm, enableDiffuseTerm, enableSpecularTerm;
};

inline float Lighting::GetLightIntensity(int lightID)
{
  return lightIntensity[lightID];
}

inline void Lighting::GetLightPosition(int lightID, float * pos)
{
  memcpy(pos, &(lightPosition[4*lightID]), 4* sizeof(float));
}

inline void Lighting::GetLightAmbientTerm(int lightID, float * Ka)
{
  memcpy(Ka, &(lKa[4*lightID]), 4* sizeof(float));
}

inline void Lighting::GetLightDiffuseTerm(int lightID, float * Kd)
{
  memcpy(Kd, &(lKd[4*lightID]), 4* sizeof(float));
}

inline void Lighting::GetLightSpecularTerm(int lightID, float * Ks)
{
  memcpy(Ks, &(lKs[4*lightID]), 4* sizeof(float));
}

inline float Lighting::GetAmbientIntensity()
{
  return ambientIntensity;
}

inline bool Lighting::IsLightEnabled(int lightID)
{
  return lightEnabled[lightID];
}

#endif

