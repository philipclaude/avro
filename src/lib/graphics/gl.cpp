#include "graphics/gl.h"

#include <stdio.h>

namespace luna
{

namespace graphics
{

int
checkOpenGLError(const char* file, int line)
{
  // Returns 1 if an OpenGL error occurred, 0 otherwise.
  GLenum glerr;
  int    retcode = 0;

  glerr = glGetError();
  while (glerr != GL_NO_ERROR)
  {
    const char* message = "";
    switch (glerr)
    {
      case GL_INVALID_ENUM:
      	message = "invalid enum";
      	break;
      case GL_INVALID_VALUE:
      	message = "invalid value";
      	break;
      case GL_INVALID_OPERATION:
      	message = "invalid operation";
      	break;
      case GL_INVALID_FRAMEBUFFER_OPERATION:
      	message = "invalid framebuffer operation";
      	break;
      case GL_OUT_OF_MEMORY:
      	message = "out of memory";
      	break;
      default:
      	message = "unknown error";
    }

    printf("glError in file %s at line %d: %s\n", file, line, message);
    retcode = 1;
    glerr = glGetError();
  }
  return retcode;
}

void
dumpGLInfo(bool dumpExtensions)
{
  const GLubyte *renderer = glGetString( GL_RENDERER );
  const GLubyte *vendor = glGetString( GL_VENDOR );
  const GLubyte *version = glGetString( GL_VERSION );
  const GLubyte *glslVersion = glGetString( GL_SHADING_LANGUAGE_VERSION );

  GLint major, minor, samples, sampleBuffers;
  glGetIntegerv(GL_MAJOR_VERSION, &major);
  glGetIntegerv(GL_MINOR_VERSION, &minor);
  glGetIntegerv(GL_SAMPLES, &samples);
  glGetIntegerv(GL_SAMPLE_BUFFERS, &sampleBuffers);

	printf("-------------------------------------------------------------\n");
  printf("GL Vendor    : %s\n", vendor);
  printf("GL Renderer  : %s\n", renderer);
  printf("GL Version   : %s\n", version);
  printf("GL Version   : %d.%d\n", major, minor);
  printf("GLSL Version : %s\n", glslVersion);
	printf("MSAA samples : %d\n", samples);
	printf("MSAA buffers : %d\n", sampleBuffers);
  printf("-------------------------------------------------------------\n");

  if (dumpExtensions)
  {
    GLint next;
    glGetIntegerv(GL_NUM_EXTENSIONS, &next);
    for (int i=0;i<next;i++)
      printf("%s\n", glGetStringi(GL_EXTENSIONS, i));
  }
}

} // graphics

} // luna
