#ifndef URSA_LIB_GRAPHICS_GL_H_
#define URSA_LIB_GRAPHICS_GL_H_

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

int checkForOpenGLError(const char*, int);
void dumpGLInfo(bool dumpExtensions = false);

} // graphics

} // ursa

#endif
