/*
  Virtual buffer object OpenGL extension function prototypes (under Windows)
*/

#ifndef _VBO_H_
#define _VBO_H_

#ifdef WIN32

//#include <windows.h>
#include <GL/gl.h>
#include <GL/glext.h>

// call this to initialize the VBO extension
// returns "true" on success and "false" on failure
// if successful, you can then use "glBindBufferARB", 
//   "glBufferDataARB" and "glGenBuffersARB", and "glDeleteBuffersARB"
bool InitializeVBOs(void);

// === no user-callable code below ===

// external function pointers 
#define NODEFINEEXTERNBUFFER
#ifndef NODEFINEEXTERNBUFFER
extern PFNGLBINDBUFFERARBPROC glBindBufferARB;
extern PFNGLBUFFERDATAARBPROC glBufferDataARB;
extern PFNGLGENBUFFERSARBPROC glGenBuffersARB;
extern PFNGLDELETEBUFFERSARBPROC glDeleteBuffersARB;
#endif

#else
  bool InitializeVBOs(void);
#endif

#endif 

