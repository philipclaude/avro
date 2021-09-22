//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_GL_H_
#define avro_LIB_GRAPHICS_GL_H_

#include "common/error.h"
#include "common/tools.h"

#include "avro_config.h"
#include "avro_types.h"

#if AVRO_WITH_GL

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

#include <map>
#include <set>
#include <vector>

namespace avro
{

namespace graphics
{

class Primitive;
class ShaderProgram;
class ShaderLibrary;

typedef GLuint  gl_index;
typedef GLfloat gl_float;

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
  avro_assert(!error); \
}

} // graphics

} // avro

#else // AVRO_WITH_GL

#include "common/error.h"

namespace avro
{

namespace graphics
{

class Primitive;
class ShaderLibrary;

typedef unsigned int gl_index;
typedef float   gl_float;
typedef gl_index GLint;
typedef gl_index GLuint;

#define GL_ELEMENT_ARRAY_BUFFER 0
#define GL_ARRAY_BUFFER 1
#define GL_PATCH_VERTICES 2
#define GL_UNSIGNED_INT 3
#define GL_LINES 4
#define GL_TRIANGLES 5
#define GL_POINTS 6
#define GL_PATCHES 11
#define GL_TEXTURE0 7
#define GL_TEXTURE_BUFFER 8
#define GL_R32F 9
#define GL_FLOAT 10

#define GLFW_PRESS 0
#define GLFW_MOD_SHIFT 1
#define GLFW_MOD_CONTROL 2
#define GLFW_RELEASE 3
#define GLFW_KEY_P 4

#define GL_CALL(X)

struct GLFWwindow {};
inline void glfwInit() { avro_assert_not_reached; }
inline void glfwWindowHint( int ) { avro_assert_not_reached; }
inline void glfwSwapBuffers( struct GLFWwindow* ) { avro_assert_not_reached; }
inline void glfwTerminate() { avro_assert_not_reached; }
inline void glfwDestroyWindow(struct GLFWwindow* ) { printf("this should not be reached\n"); }
inline void glfwGetCursorPos( struct GLFWwindow* , double* , double* ) { avro_assert_not_reached; }

inline void glGenBuffers( gl_index , gl_index* ) { avro_assert_not_reached; }
inline void glDeleteBuffers( gl_index , gl_index* ) { printf("this should not be reached\n"); }
inline void glDeleteTextures( gl_index , gl_index* ) { printf("this should not be reached\n"); }
inline void glDeleteVertexArrays( gl_index , gl_index* ) { printf("this should not be reached\n"); }
inline void glActiveTexture( gl_index ) { avro_assert_not_reached; }
inline void glUniform1i( gl_index , int ) { avro_assert_not_reached; }
inline gl_index glGetUniformLocation( gl_index , const std::string& ) { avro_assert_not_reached; }
inline void glClearColor( float , float , float , float ) { avro_assert_not_reached; }
inline void glViewport( int , int , int , int ) { avro_assert_not_reached; }

#define glEnable(x)
#define glDisable(x)
#define glClear(x)

class ShaderProgram {
public:
  ShaderProgram( const std::string& , bool , const std::vector<std::string> );
  void use() {}
  template<typename T> void setUniform( const std::string& , const T& ) { avro_assert_not_reached; }

  gl_index handle() const { return 0; }
  bool has_tessellation_shader() const { return false; }
};

}

}

#endif // AVRO_WITH_GL

#endif
