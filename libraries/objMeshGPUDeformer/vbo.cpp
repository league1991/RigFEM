#include "stdafx.h"
#include "vbo.h"

#ifdef WIN32

// PFNGLBINDBUFFERARBPROC glBindBufferARB = NULL;
// PFNGLBUFFERDATAARBPROC glBufferDataARB = NULL;
// PFNGLGENBUFFERSARBPROC glGenBuffersARB = NULL;
// PFNGLDELETEBUFFERSARBPROC glDeleteBuffersARB = NULL;

bool InitializeVBOs(void)
{
  if (glBindBufferARB == NULL)
    glBindBufferARB = (PFNGLBINDBUFFERARBPROC) wglGetProcAddress("glBindBufferARB");

  if (glBufferDataARB == NULL)
    glBufferDataARB = (PFNGLBUFFERDATAARBPROC) wglGetProcAddress("glBufferDataARB");

  if (glGenBuffersARB == NULL)
    glGenBuffersARB = (PFNGLGENBUFFERSARBPROC) wglGetProcAddress("glGenBuffersARB");

  if (glDeleteBuffersARB == NULL)
    glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) wglGetProcAddress("glDeleteBuffersARB");

  return (glBindBufferARB && glBufferDataARB && glGenBuffersARB && glDeleteBuffersARB);
}

#else
  bool InitializeVBOs(void) { return true; }
#endif

