//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/gl.h"

#include <stdio.h>

namespace avro
{

namespace graphics
{

#if AVRO_WITH_GL

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
dumpGLInfo(bool dumpExtensions) {

  const GLubyte *renderer = glGetString( GL_RENDERER );
  const GLubyte *vendor = glGetString( GL_VENDOR );
  const GLubyte *version = glGetString( GL_VERSION );
  const GLubyte *glslVersion = glGetString( GL_SHADING_LANGUAGE_VERSION );

  GLint major, minor;
  glGetIntegerv(GL_MAJOR_VERSION, &major);
  glGetIntegerv(GL_MINOR_VERSION, &minor);

	printf("-------------------------------------------------------------\n");
  printf("--> vendor: %s\n", vendor);
  printf("--> device: %s\n", renderer);
  printf("--> driver: %s (@ OpenGL %d.%d)\n", version, major, minor);
  printf("--> @ glsl: %s\n", glslVersion);
  printf("-------------------------------------------------------------\n");

  if (dumpExtensions) {
    GLint nb_ext;
    glGetIntegerv(GL_NUM_EXTENSIONS, &nb_ext);
    for (int i = 0; i < nb_ext; i++)
      printf("%s\n", glGetStringi(GL_EXTENSIONS, i));
  }
}

#endif

} // graphics

} // avro
