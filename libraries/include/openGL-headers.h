#if defined(WIN32) || defined(linux)
  #include <GL/gl.h> 
  #include <GL/glu.h> 
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
  #include <GLUT/glut.h>
#endif

