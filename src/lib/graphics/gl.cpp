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
#include "graphics/primitive.h"
#include "graphics/scene.h"
#include "graphics/shader.h"

#include <stdio.h>

namespace avro
{

namespace graphics
{

#ifdef AVRO_WITH_GL

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

OpenGL_Manager::OpenGL_Manager()
{
  shaders_ = std::make_shared<ShaderLibrary>();
}

void
OpenGL_Manager::select_shader( Primitive& primitive , const std::string& name )
{
  shader_[&primitive] = &(*shaders_)[name];
}

void
OpenGL_Manager::select_shader( const std::string& name , const std::string& shader_name )
{
  aux_vao_shader_[name] = &(*shaders_)[shader_name];
}

void
OpenGL_Manager::create_shaders()
{
  shaders_->create();
}

void
OpenGL_Manager::write( Primitive& primitive )
{

  //if (vao_edges_.find(&primitive)!=vao_edges_.end()) return;

  GLuint vbo_position;
  GLuint vbo_colour;
  GLuint vbo_normal;
  GLuint vbo_triangles;
  GLuint vbo_edges;

  GLuint vbo_feedback_points;
  GLuint vbo_feedback_edges;
  GLuint vbo_feedback_triangles;

  if (primitive_.find(&primitive)==primitive_.end())
  {
    // generate new buffers
    std::vector<GLuint> vbo(9);
    GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );

    vbo_position  = vbo[0];
    vbo_colour    = vbo[1];
    vbo_normal    = vbo[2];
    vbo_triangles = vbo[3];
    vbo_edges     = vbo[4];

    vbo_feedback_points    = vbo[6];
    vbo_feedback_edges     = vbo[7];
    vbo_feedback_triangles = vbo[8];

    // keep track of the vbo's in case we need to update them
    vbo_points_.insert( {&primitive,vbo_position} );
    vbo_edges_.insert( {&primitive,vbo_edges} );
    vbo_triangles_.insert( {&primitive,vbo_triangles} );
    vbo_normals_.insert( {&primitive,vbo_normal} );
    vbo_colors_.insert( {&primitive,vbo_colour} );

    vbo_feedback_points_.insert( {&primitive,vbo_feedback_points} );
    vbo_feedback_edges_.insert( {&primitive,vbo_feedback_edges} );
    vbo_feedback_triangles_.insert( {&primitive,vbo_feedback_triangles} );

    primitive_.insert( &primitive );
  }
  else
  {
    // retrieve the buffer indices, so we can overwrite them
    avro_assert( vbo_points_.find(&primitive) != vbo_points_.end() );
    avro_assert( vbo_edges_.find(&primitive) != vbo_edges_.end() );
    avro_assert( vbo_triangles_.find(&primitive) != vbo_triangles_.end() );
    avro_assert( vbo_normals_.find(&primitive) != vbo_normals_.end() );
    avro_assert( vbo_colors_.find(&primitive) != vbo_colors_.end() );
    avro_assert( vbo_feedback_points_.find(&primitive) != vbo_feedback_points_.end() );
    avro_assert( vbo_feedback_edges_.find(&primitive) != vbo_feedback_edges_.end() );
    avro_assert( vbo_feedback_triangles_.find(&primitive) != vbo_feedback_triangles_.end() );

    vbo_position  = vbo_points_[&primitive];
    vbo_edges     = vbo_edges_[&primitive];
    vbo_triangles = vbo_triangles_[&primitive];
    vbo_normal    = vbo_normals_[&primitive];
    vbo_colour    = vbo_colors_[&primitive];

    vbo_feedback_points    = vbo_feedback_points_[&primitive];
    vbo_feedback_edges     = vbo_feedback_edges_[&primitive];
    vbo_feedback_triangles = vbo_feedback_triangles_[&primitive];
  }


  if (primitive.number()>=2)
  {
    // bind the triangles
    std::vector<GLuint> triangles( primitive.triangles().begin() , primitive.triangles().end() );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_triangles) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(GLuint), triangles.data() , GL_STATIC_DRAW) );
  }

  // bind the position buffer
  std::vector<GLfloat> points( primitive.points().begin() , primitive.points().end() );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat), points.data() , GL_STATIC_DRAW) );

  // bind the transform feedback buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_feedback_points) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat) , NULL , GL_STATIC_COPY) );
  glBindBuffer(GL_ARRAY_BUFFER,0);

  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_feedback_edges) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat) , NULL , GL_STATIC_COPY) );
  glBindBuffer(GL_ARRAY_BUFFER,0);

  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_feedback_triangles) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat) , NULL , GL_STATIC_COPY) );
  glBindBuffer(GL_ARRAY_BUFFER,0);

  // bind the colour buffer
  std::vector<GLfloat> colors( primitive.colors().begin() , primitive.colors().end() );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_colour) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(GLfloat), colors.data() , GL_STATIC_DRAW) );

  // bind the normal buffer
  if (primitive.number()>=2)
  {
    std::vector<GLfloat> normals( primitive.normals().begin() , primitive.normals().end() );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_normal) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(GLfloat), normals.data() , GL_STATIC_DRAW) );
  }

  // bind the triangle data
  GLuint id_vao_triangles;
  GL_CALL( glGenVertexArrays( 1, &id_vao_triangles ) );
  GL_CALL( glBindVertexArray(id_vao_triangles) );
  vao_triangles_.insert( {&primitive,id_vao_triangles} );

  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_triangles ) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_colour ) );
  GL_CALL( glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(1) );

  if (primitive.number()>=2)
  {
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_normal ) );
    GL_CALL( glVertexAttribPointer( 2, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(2) );
  }

  // bind the edges
  GLuint id_vao_edges;
  GL_CALL( glGenVertexArrays( 1, &id_vao_edges ) );
  GL_CALL( glBindVertexArray(id_vao_edges) );
  vao_edges_.insert( {&primitive,id_vao_edges} );

  std::vector<GLuint> edges( primitive.edges().begin() , primitive.edges().end() );
  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_edges ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges.size() * sizeof(GLuint), edges.data() , GL_STATIC_DRAW) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // bind the points
  GLuint id_vao_points;
  GL_CALL( glGenVertexArrays( 1, &id_vao_points ) );
  GL_CALL( glBindVertexArray(id_vao_points) );
  vao_points_.insert( {&primitive,id_vao_points} );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // bind the transform feedback buffer
  GLuint id_vao_feedback_points;
  GL_CALL( glGenVertexArrays( 1, &id_vao_feedback_points ) );
  GL_CALL( glBindVertexArray( id_vao_feedback_points ) );
  vao_feedback_points_.insert( {&primitive,id_vao_feedback_points} );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , vbo_feedback_points ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GLuint id_vao_feedback_edges;
  GL_CALL( glGenVertexArrays( 1, &id_vao_feedback_edges ) );
  GL_CALL( glBindVertexArray( id_vao_feedback_edges ) );
  vao_feedback_edges_.insert( {&primitive,id_vao_feedback_edges} );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , vbo_feedback_edges ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GLuint id_vao_feedback_triangles;
  GL_CALL( glGenVertexArrays( 1, &id_vao_feedback_triangles ) );
  GL_CALL( glBindVertexArray( id_vao_feedback_triangles) );
  vao_feedback_triangles_.insert( {&primitive,id_vao_feedback_triangles} );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , vbo_feedback_triangles ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

void
OpenGL_Manager::write( const std::string& name , coord_t number , const std::vector<real_t>& points0 , const std::vector<index_t>& edges0 , const std::vector<index_t>& triangles0 , const std::vector<real_t>& colors0 )
{
  GLuint vbo_points;
  GLuint vbo_color;
  GLuint vbo_triangles;
  GLuint vbo_edges;

  // generate new buffers
  std::vector<GLuint> vbo(4);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );

  vbo_points    = vbo[0];
  vbo_color     = vbo[1];
  vbo_triangles = vbo[2];
  vbo_edges     = vbo[3];

  // write the points
  std::vector<GLfloat> points( points0.begin() , points0.end() );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_points) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(GLfloat), points.data() , GL_STATIC_DRAW) );

  // write to the triangles buffer
  if (number>=2)
  {
    // write the triangle data
    std::vector<GLuint> triangles( triangles0.begin() , triangles0.end() );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_triangles) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(GLuint), triangles.data() , GL_STATIC_DRAW) );

    // write the color data
    std::vector<GLfloat> colors( colors0.begin() , colors0.end() );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, vbo_color) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(GLfloat), colors.data() , GL_STATIC_DRAW) );

    // bind the triangle data
    GLuint id_vao_triangles;
    GL_CALL( glGenVertexArrays( 1, &id_vao_triangles ) );
    GL_CALL( glBindVertexArray(id_vao_triangles) );

    GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_triangles ) );

    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_points ) );
    GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(0) );

    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_color ) );
    GL_CALL( glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(1) );

    std::string id_name = name + "-triangles";
    if (aux_vao_.find(id_name)==aux_vao_.end())
    {
      aux_vao_.insert( {id_name,id_vao_triangles} );
      aux_vao_size_.insert( {id_name,triangles.size() } );
    }
    else
    {
      aux_vao_[id_name] = id_vao_triangles;
      aux_vao_size_[id_name] = triangles.size();
    }
  }

  // write to the edges buffer
  if (number>=1)
  {
    GLuint id_vao_edges;
    GL_CALL( glGenVertexArrays( 1, &id_vao_edges ) );
    GL_CALL( glBindVertexArray(id_vao_edges) );

    std::vector<GLuint> edges( edges0.begin() , edges0.end() );
    GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, vbo_edges ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges.size() * sizeof(GLuint), edges.data() , GL_STATIC_DRAW) );

    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_points ) );
    GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(0) );

    // bind the points
    GLuint id_vao_points;
    GL_CALL( glGenVertexArrays( 1, &id_vao_points ) );
    GL_CALL( glBindVertexArray(id_vao_points) );

    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_points ) );
    GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(0) );

    if (number>=1)
    {
      std::string id_name = name + "-edges";
      if (aux_vao_.find(id_name)==aux_vao_.end())
      {
        aux_vao_.insert( {id_name,id_vao_edges} );
        aux_vao_size_.insert( {id_name,edges.size() } );
      }
      else
      {
        aux_vao_[id_name] = id_vao_edges;
        aux_vao_size_[id_name] = edges.size();
      }
    }
  }

  // write to the points buffer
  GLuint id_vao_points;
  GL_CALL( glGenVertexArrays( 1, &id_vao_points ) );
  GL_CALL( glBindVertexArray(id_vao_points) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, vbo_points ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  select_shader(name , "wv" );

}

void
OpenGL_Manager::draw( Primitive& primitive , TransformFeedbackResult* feedback )
{
  // indicate to the gl that we want to use the shader
  avro_assert_msg( shader_.find(&primitive)!=shader_.end() , "shader not set" );
  ShaderProgram shader = *shader_.at(&primitive);
  shader.use();

  GLuint query;
  GLuint buffer;
  if (feedback!=nullptr)
  {

    glGenQueries(1, &query );
    GL_CALL( glBeginQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN , query ) );

    GL_CALL( glGenBuffers( 1 , &buffer ) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glBufferData( GL_TRANSFORM_FEEDBACK_BUFFER , 24*primitive.triangles().size()*sizeof(GLfloat) , NULL , GL_DYNAMIC_COPY ) );
    GL_CALL( glBindBufferBase( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , buffer ) );
    GL_CALL( glBindVertexArray( vao_feedback_triangles_[&primitive] ) );

    GL_CALL( glEnable(GL_RASTERIZER_DISCARD) );
    GL_CALL( glBeginTransformFeedback(GL_TRIANGLES) );
  }

  // bind the vao associated with this primitive
  if (primitive.number()>=2)
  {
    // set the transparency
    shader.setUniform( "xpar" , (float)primitive.transparency() );

    // draw the triangles
    if (primitive.triangles_on())
    {
      GL_CALL( glBindVertexArray( vao_triangles_.at(&primitive) ) );
      GL_CALL( glDrawElements(GL_TRIANGLES, primitive.triangles().size() , GL_UNSIGNED_INT , 0 ) )
    }
  }

  if (feedback!=nullptr)
  {
    GL_CALL( glEndTransformFeedback() );
    GL_CALL( glDisable(GL_RASTERIZER_DISCARD) );

    GLuint nb_primitives = 0;
    GL_CALL( glEndQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN ) );
    GL_CALL( glGetQueryObjectuiv( query , GL_QUERY_RESULT , &nb_primitives ) );

    // there are 4 coordinates for every output and 3 vertices per prim (color + coord)
    index_t nb_data_per_prim = 4*3*2;
    index_t size = nb_primitives*nb_data_per_prim;

    GLfloat* feedbackBuffer = (GLfloat*) malloc( size*sizeof(GLfloat) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glGetBufferSubData( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , size*sizeof(GLfloat) , feedbackBuffer ) );

    feedback->append_triangles(feedbackBuffer,nb_primitives);

    // begin transform feedback for edges
    GL_CALL( glBeginQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN , query ) );

    GL_CALL( glGenBuffers( 1 , &buffer ) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glBufferData( GL_TRANSFORM_FEEDBACK_BUFFER , 16*primitive.edges().size()*sizeof(GLfloat) , NULL , GL_DYNAMIC_COPY ) );
    GL_CALL( glBindBufferBase( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , buffer ) );
    GL_CALL( glBindVertexArray( vao_feedback_edges_[&primitive] ) );

    GL_CALL( glEnable(GL_RASTERIZER_DISCARD) );
    GL_CALL( glBeginTransformFeedback(GL_LINES) );
  }

  if (primitive.number()>=1)
  {
    if (primitive.edges_on())
    {
      GL_CALL( glBindVertexArray( vao_edges_.at(&primitive) ) );
      GL_CALL( glDrawElements( GL_LINES , primitive.edges().size() , GL_UNSIGNED_INT , 0 ) )
    }
  }

  if (feedback!=nullptr)
  {
    GL_CALL( glEndTransformFeedback() );
    GL_CALL( glDisable(GL_RASTERIZER_DISCARD) );

    GLuint nb_primitives = 0;
    GL_CALL( glEndQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN ) );
    GL_CALL( glGetQueryObjectuiv( query , GL_QUERY_RESULT , &nb_primitives ) );

    // there are 4 coordinates for every output and 3 vertices per prim (color + coord)
    index_t nb_data_per_prim = 4*2*2;
    index_t size = nb_primitives*nb_data_per_prim;

    GLfloat* feedbackBuffer = (GLfloat*) malloc( size*sizeof(GLfloat) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glGetBufferSubData( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , size*sizeof(GLfloat) , feedbackBuffer ) );

    feedback->append_edges(feedbackBuffer,nb_primitives);
  }


  // draw the points
  if (primitive.points_on())
  {
    GL_CALL( glBindVertexArray(vao_points_.at(&primitive) ) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , primitive.points().size()/3 ) );
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );

  for (index_t k=0;k<primitive.nb_children();k++)
  {
    draw( primitive.child(k) , feedback );
  }
}

void
OpenGL_Manager::draw( SceneGraph& scene , TransformFeedbackResult* feedback )
{
  if (scene.update())
    shaders_->set_matrices(scene);

  // draw the primitives
  for (index_t k=0;k<scene.nb_primitives();k++)
  {
    draw( scene.primitive(k) , feedback );
  }

  scene.set_update(false);
}

void
OpenGL_Manager::draw( const std::string& name , coord_t number , const DrawingParameters& params )
{

  std::string id_triangles = name + "-triangles";
  std::string id_edges = name + "-edges";

  ShaderProgram* shader = aux_vao_shader_[name];
  shader->use();
  shader->setUniform( "xpar" , params.transparency );
  shader->setUniform( "MVP" , params.mvp );
  shader->setUniform( "wLight" , params.lighting );

  if (number>=1)
  {
    std::map<std::string,GLuint>::iterator it;
    it = aux_vao_.find(id_edges);
    avro_assert_msg( it!=aux_vao_.end() , "could not find vao associated with %s" , id_edges.c_str() );

    GLuint vao = it->second;
    index_t size = aux_vao_size_[id_edges];

    GL_CALL( glBindVertexArray( vao ) );
    GL_CALL( glDrawElements( GL_LINES ,size , GL_UNSIGNED_INT , 0 ) );
  }
  if (number==2)
  {
    std::map<std::string,GLuint>::iterator it;
    it = aux_vao_.find(id_triangles);
    avro_assert_msg( it!=aux_vao_.end() , "could not find vao associated with %s" , id_triangles.c_str() );

    GLuint vao = it->second;
    index_t size = aux_vao_size_[id_triangles];

    GL_CALL( glBindVertexArray( vao ) );
    GL_CALL( glDrawElements( GL_TRIANGLES ,size , GL_UNSIGNED_INT , 0 ) );
  }

  GL_CALL( glBindVertexArray(0) );
}

#endif

} // graphics

} // avro
