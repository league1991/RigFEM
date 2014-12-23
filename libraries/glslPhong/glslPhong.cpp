/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "glslPhong" library , Copyright (C) 2014 USC                          *
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

/*
   Phong shading ("per-pixel lighting") in GLSL.
   Requires OpenGL 2.0 or higher
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(WIN32) || defined(linux)
  #define USE_GLEW
#elif defined(__APPLE__)
  #include "openGL-headers.h"
  #include <OpenGL/glext.h>
#endif

#ifdef USE_GLEW
  #include <GL/glew.h>
#endif

#include "glslPhong.h"

char GLSLPhong::vertexShaderStringAll [] = "varying vec3 normal, eyeVec;\n"
"#define MAX_LIGHTS 8\n"
"varying vec3 lightDir[MAX_LIGHTS];\n"
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"  normal = gl_NormalMatrix * gl_Normal;\n"
"  vec4 vVertex = gl_ModelViewMatrix * gl_Vertex;\n"
"  eyeVec = -vVertex.xyz;\n"
"  int i;\n"
"  for (i=0; i<MAX_LIGHTS; ++i)\n"
"    lightDir[i] = vec3(gl_LightSource[i].position.xyz - vVertex.xyz);\n"
"}\n";

char GLSLPhong::fragmentShaderStringAll [] = "varying vec3 normal, eyeVec;\n"
"#define MAX_LIGHTS 8\n"
"varying vec3 lightDir[MAX_LIGHTS];\n"
"void main (void)\n"
"{\n"
"  vec4 final_color = gl_FrontLightModelProduct.sceneColor;\n"
"  vec3 N = normalize(normal);\n"
"  int i;\n"
"  for (i=0; i<MAX_LIGHTS; ++i)\n"
"  {\n"
"    final_color += gl_LightSource[i].ambient * gl_FrontMaterial.ambient;\n"
"    vec3 L = normalize(lightDir[i]);\n"
"    float lambertTerm = dot(N,L);\n"
"    if (lambertTerm > 0.0)\n"
"    {\n"
"      final_color += gl_LightSource[i].diffuse * gl_FrontMaterial.diffuse * lambertTerm;\n"
"      vec3 E = normalize(eyeVec);\n"
"      vec3 R = reflect(-L, N);\n"
"      float specular = pow(max(dot(R, E), 0.0), gl_FrontMaterial.shininess);\n"
"      final_color += gl_LightSource[i].specular * gl_FrontMaterial.specular * specular;\n"
"    }\n"
"  }\n"
"  gl_FragColor = final_color;\n"
"}\n";

//"  gl_Position = ftransform();\n"
//"  gl_Position = projection_matrix * modelview_matrix * vec4(vertex, 1.0);\n"
//"uniform Transformation {\n"
//"    mat4 projection_matrix;\n"
//"    mat4 modelview_matrix;\n"
//"};\n"
//"  lightDir[%d] = vec3(gl_LightSource[%d].position.xyz - vVertex.xyz);\n";

char GLSLPhong::vertexShaderStringPrologue [] = ""
"varying vec3 normal, eyeVec;\n"
"#define MAX_LIGHTS 8\n"
"varying vec3 lightDir[MAX_LIGHTS];\n"
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"  normal = gl_NormalMatrix * gl_Normal;\n"
"  vec4 vVertex = gl_ModelViewMatrix * gl_Vertex;\n"
"  eyeVec = -vVertex.xyz;\n"
"  vec4 transformedLightPos;\n";

char GLSLPhong::vertexShaderStringCore [] = ""
"  transformedLightPos = vec4(gl_LightSource[%d].position.xyz, 1);\n"
"  lightDir[%d] = vec3(transformedLightPos.xyz - vVertex.xyz);\n";

char GLSLPhong::vertexShaderStringEpilogue [] = "}\n";

char GLSLPhong::fragmentShaderStringPrologue [] = ""
"varying vec3 normal, eyeVec;\n"
"#define MAX_LIGHTS 8\n"
"varying vec3 lightDir[MAX_LIGHTS];\n"
"void main (void)\n"
"{\n"
"  vec4 final_color = gl_FrontLightModelProduct.sceneColor;\n"
"  vec3 N = normalize(normal);\n";

/*
// no backspace culling
char GLSLPhong::fragmentShaderStringCore [] = ""
"  {\n"
"    final_color += gl_LightSource[%d].ambient * gl_FrontMaterial.ambient;\n"
"    vec3 L = normalize(lightDir[%d]);\n"
"    float lambertTerm = dot(N,L);\n"
"    if (lambertTerm > 0.0)\n"
"    {\n"
"      final_color += gl_LightSource[%d].diffuse * gl_FrontMaterial.diffuse * lambertTerm;\n"
"      vec3 E = normalize(eyeVec);\n"
"      vec3 R = reflect(-L, N);\n"
"      float specular = pow(max(dot(R, E), 0.0), gl_FrontMaterial.shininess);\n"
"      final_color += gl_LightSource[%d].specular * gl_FrontMaterial.specular * specular;\n"
"    }\n"
"  }\n";
*/

// with backspace culling
char GLSLPhong::fragmentShaderStringCore [] = ""
"  {\n"
"    final_color += gl_LightSource[%d].ambient * gl_FrontMaterial.ambient;\n"
"    vec3 L = normalize(lightDir[%d]);\n"
"    float lambertTerm = dot(N,L);\n"
"    float viewingTerm = dot(N,eyeVec);\n"
"    if ((lambertTerm > 0.0) && (viewingTerm > 0.0))\n"
"    {\n"
"      final_color += gl_LightSource[%d].diffuse * gl_FrontMaterial.diffuse * lambertTerm;\n"
"      vec3 E = normalize(eyeVec);\n"
"      vec3 R = reflect(-L, N);\n"
"      float specular = pow(max(dot(R, E), 0.0), gl_FrontMaterial.shininess);\n"
"      final_color += gl_LightSource[%d].specular * gl_FrontMaterial.specular * specular;\n"
"    }\n"
"  }\n";

// no backspace culling
//"    if (lambertTerm > 0.0)\n"

// with backface culling:
//"    float viewingTerm = dot(N,eyeVec);\n"
//"    if ((lambertTerm > 0.0) && (viewingTerm > 0.0))\n"

char GLSLPhong::fragmentShaderStringEpilogue [] = ""
"  gl_FragColor = final_color;\n"
"}\n";

void GLSLPhong::checkError(GLint status, const char *msg)
{
  if (status != GL_TRUE)
  {
    printf("%s\n", msg);
    exit(1);
  }
}

void GLSLPhong::checkShaderCompilation(GLint shaderID)
{
  GLint blen = 0; 
  GLsizei slen = 0;

  glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH , &blen);       

  if (blen > 1)
  {
    GLchar * compiler_log = (GLchar*) malloc (blen);
    glGetShaderInfoLog(shaderID, blen, &slen, compiler_log);
    printf("compiler_log:\n%s\n", compiler_log);
    free(compiler_log);
  }
}

GLSLPhong::GLSLPhong()
{
  #ifdef USE_GLEW
    GLenum err = glewInit();
    if (GLEW_OK != err)
    {
      /* Problem: glewInit failed, something is seriously wrong. */
      fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
  #endif

  for(int programID=0; programID<256; programID++)
  {
    //printf("Preparing program %d...\n", programID);
    programs[programID] = glCreateProgram();
    //printf("Created program %d.\n", programID);

    char vertexShaderString[20000];
    char fragmentShaderString[20000];
    char buffer[20000];

    strcpy(vertexShaderString, vertexShaderStringPrologue);
    strcpy(fragmentShaderString, fragmentShaderStringPrologue);
    for(int light=0; light<8; light++)
    {
      if (((programID >> light) & 1) == 0)
        continue;
      sprintf(buffer, vertexShaderStringCore, light, light);
      strcat(vertexShaderString, buffer);
      sprintf(buffer, fragmentShaderStringCore, light, light, light, light);
      strcat(fragmentShaderString, buffer);
    }
    strcat(vertexShaderString, vertexShaderStringEpilogue);
    strcat(fragmentShaderString, fragmentShaderStringEpilogue);

    const GLchar * vertexSource = vertexShaderString;
    const GLchar * fragmentSource = fragmentShaderString;

    //printf("programID=%d\n", programID);
    //printf("vertexSource:\n%s\n\n", vertexSource);
    //printf("fragmentSource:\n%s\n", fragmentSource);

    //if (programID == 128)
      //exit(1);

    /* create program and shader objects */

    vertexShaders[programID] = glCreateShader(GL_VERTEX_SHADER);
    fragmentShaders[programID] = glCreateShader(GL_FRAGMENT_SHADER);

    /* attach shaders to the program object */

    glAttachShader(programs[programID], vertexShaders[programID]);
    glAttachShader(programs[programID], fragmentShaders[programID]);

    /* load shaders */

    glShaderSource(vertexShaders[programID], 1, &vertexSource, NULL);
    glShaderSource(fragmentShaders[programID], 1, &fragmentSource, NULL);

    /* compile shaders */

    glCompileShader(vertexShaders[programID]);
    glCompileShader(fragmentShaders[programID]);

    /* error check */

    GLint status = 0;
    glGetShaderiv(vertexShaders[programID], GL_COMPILE_STATUS, &status);
    checkShaderCompilation(vertexShaders[programID]);
    checkError(status, "Failed to compile the vertex shader.");

    glGetShaderiv(fragmentShaders[programID], GL_COMPILE_STATUS, &status);
    checkShaderCompilation(fragmentShaders[programID]);
    checkError(status, "Failed to compile the fragment shader.");

    /* set up uniform parameter */
    //timeParam = glGetUniformLocation(program, "time");

    /* link */

    glLinkProgram(programs[programID]);

    glGetShaderiv(programs[programID], GL_LINK_STATUS, &status);
    checkError(status, "Failed to link the shader program object.");
  }

  Disable();
}

GLSLPhong::~GLSLPhong()
{
  for(int programID=0; programID<256; programID++)
  {
    glDetachShader(programs[programID], fragmentShaders[programID]);
    glDetachShader(programs[programID], vertexShaders[programID]);
    glDeleteShader(fragmentShaders[programID]);
    glDeleteShader(vertexShaders[programID]);
    glDeleteProgram(programs[programID]);
  }
}

void GLSLPhong::Enable()
{
  GLenum lightArray[8] = { 
    GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3,
    GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 };
  int programID = 0;
  for(int light=7; light>=0; light--)
  {
    programID *= 2;
    if (glIsEnabled(lightArray[light]) == GL_TRUE)
      programID += 1;
  }
  //printf("programID=%d\n", programID);
  glUseProgram(programs[programID]);
}

void GLSLPhong::Disable()
{
  glUseProgram(0);
}

