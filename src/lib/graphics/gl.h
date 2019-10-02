#ifndef URSA_LIB_GRAPHICS_GL_H_
#define URSA_LIB_GRAPHICS_GL_H_

#include "common/error.h"

#include <glad/include/glad/glad.h>
#define GLFW_INCLUDE_NONE

#ifdef __APPLE__
#define __gl_h_
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#include <GLFW/glfw3.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#endif

namespace ursa
{

namespace graphics
{

int checkOpenGLError(const char*, int);
void dumpGLInfo(bool dumpExtensions = false);

#define GL_CALL( X ) \
{  \
  (X); \
  GLenum glerr; bool error = false; \
  glerr = glGetError(); \
  while (glerr != GL_NO_ERROR) \
  { \
    const char* message = ""; \
    switch (glerr) \
    { \
      case GL_INVALID_ENUM: \
      	message = "invalid enum"; \
      	break;  \
      case GL_INVALID_VALUE: \
      	message = "invalid value"; \
      	break; \
      case GL_INVALID_OPERATION: \
      	message = "invalid operation"; \
      	break; \
      case GL_INVALID_FRAMEBUFFER_OPERATION: \
      	message = "invalid framebuffer operation"; \
      	break; \
      case GL_OUT_OF_MEMORY: \
      	message = "out of memory"; \
      	break; \
      default: \
      	message = "unknown error"; \
    } \
    printf("OpenGL error in file %s at line %d: %s\n", __FILE__,__LINE__, message); \
    glerr = glGetError(); \
    error = true; \
  } \
  ursa_assert(!error); \
}


} // graphics

} // ursa

#endif
