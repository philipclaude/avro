#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/scene.h"
#include "graphics/shader.h"

#include <stdio.h>

namespace avro
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
OpenGL_Manager::create_shaders()
{
  shaders_->create();
}

void
OpenGL_Manager::write( Primitive& primitive )
{
  // allocate the buffers
  std::vector<GLuint> vbo(9);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );
  GLuint& vbo_position  = vbo[0];
  GLuint& vbo_colour    = vbo[1];
  GLuint& vbo_normal    = vbo[2];
  GLuint& vbo_triangles = vbo[3];
  GLuint& vbo_edges     = vbo[4];
  GLuint& vbo_points    = vbo[5]; UNUSED(vbo_points);

  GLuint& vbo_feedback_points    = vbo[6];
  GLuint& vbo_feedback_edges     = vbo[7];
  GLuint& vbo_feedback_triangles = vbo[8];

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
    printf("triangles written = %u\n",nb_primitives);

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
    printf("edges written = %u\n",nb_primitives);

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
    GL_CALL( glDrawArrays( GL_POINTS , 0 , primitive.points().size() ) );
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

  for (index_t k=0;k<scene.nb_primitives();k++)
  {
    draw( scene.primitive(k) , feedback );
  }

  scene.set_update(false);
}

} // graphics

} // avro
