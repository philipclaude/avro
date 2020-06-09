//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_GL_H_
#define avro_LIB_GRAPHICS_GL_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

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

#include "graphics/manager.h"
#include "graphics/shader.h"

namespace avro
{

namespace graphics
{

class Primitive;
class ShaderProgram;
class ShaderLibrary;

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

class TransformFeedbackResult
{
public:

  void append_triangles( const GLfloat* buffer , index_t nb_triangles )
  {
    index_t n = 0;
    for (index_t k=0;k<nb_triangles;k++)
    {
      //printf("primitve %lu\n",k);
      for (index_t j=0;j<3;j++)
      {
        for (coord_t d=0;d<3;d++)
          triangle_points_.push_back( buffer[n+d] );
        for (coord_t d=0;d<3;d++)
          triangle_colors_.push_back( buffer[n+4+d] );

        /*printf("v = (%g,%g,%g,%g), c = (%g,%g,%g,%g)\n",
                buffer[n  ],buffer[n+1],buffer[n+2],buffer[n+3] ,
                buffer[n+4],buffer[n+5],buffer[n+6],buffer[n+7]);*/
        n += 8;
      }
    }
  }

  void append_edges( const GLfloat* buffer , index_t nb_edges )
  {
    index_t n = 0;
    for (index_t k=0;k<nb_edges;k++)
    {
      //printf("primitve %lu\n",k);
      for (index_t j=0;j<2;j++)
      {
        for (coord_t d=0;d<3;d++)
          edge_points_.push_back( buffer[n+d] );
        for (coord_t d=0;d<3;d++)
          edge_colors_.push_back( buffer[n+4+d] );

        /*printf("v = (%g,%g,%g,%g), c = (%g,%g,%g,%g)\n",
                buffer[n  ],buffer[n+1],buffer[n+2],buffer[n+3] ,
                buffer[n+4],buffer[n+5],buffer[n+6],buffer[n+7]);*/
        n += 8;
      }
    }
  }

  const std::vector<real_t>& triangle_points() const { return triangle_points_; }
  const std::vector<real_t>& triangle_colors() const { return triangle_colors_; }

  const std::vector<real_t>& edge_points() const { return edge_points_; }
  const std::vector<real_t>& edge_colors() const { return edge_colors_; }

private:
  std::vector<real_t> triangle_points_;
  std::vector<real_t> triangle_colors_;

  std::vector<real_t> edge_points_;
  std::vector<real_t> edge_colors_;

  std::vector<real_t> vertex_points_;
  std::vector<real_t> vertex_colors_;
};

class OpenGL_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  OpenGL_Manager();

public:
  void write( Primitive& primitive );

  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr );

  void select_shader( Primitive& primitive , const std::string& name );
  void create_shaders();
  ShaderLibrary& shaders() { return *shaders_.get(); }

private:

  void set_matrices( SceneGraph& scene );
  void draw( Primitive& primitive , TransformFeedbackResult* feedback=nullptr );

  // map from primitive to vbo
  std::set<Primitive*> primitive_;

  std::map<Primitive*,index_t> vbo_points_;
  std::map<Primitive*,index_t> vbo_edges_;
  std::map<Primitive*,index_t> vbo_triangles_;
  std::map<Primitive*,index_t> vbo_colors_;
  std::map<Primitive*,index_t> vbo_normals_;

  std::map<Primitive*,index_t> vbo_feedback_points_;
  std::map<Primitive*,index_t> vbo_feedback_edges_;
  std::map<Primitive*,index_t> vbo_feedback_triangles_;

  // map from primitive to vao
  std::map<Primitive*,index_t> vao_points_;
  std::map<Primitive*,index_t> vao_edges_;
  std::map<Primitive*,index_t> vao_triangles_;

  std::map<Primitive*,index_t> vao_feedback_points_;
  std::map<Primitive*,index_t> vao_feedback_edges_;
  std::map<Primitive*,index_t> vao_feedback_triangles_;

  std::map<Primitive*,ShaderProgram*> shader_;
  std::shared_ptr<ShaderLibrary> shaders_;
};


} // graphics

} // avro

#endif
