#include "graphics/new/application.h"
#include "graphics/new/managers.h"
#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

namespace avro
{

namespace graphics
{

OpenGL4_Manager::OpenGL4_Manager() {

  // initialize GLFW
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

}

void
OpenGL4_Manager::write( VertexAttributeObject& vao ) {

  gl_index& vertex_array = vao.vertex_array();

  GL_CALL( glGenVertexArrays( 1, &vertex_array ) );
  GL_CALL( glBindVertexArray(vertex_array) );
  vao_.push_back(vertex_array);

  gl_index& points_buffer = vao.points().buffer();
  GL_CALL( glGenBuffers( 1 , &points_buffer ) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, points_buffer ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * vao.points().coordinates().size() , vao.points().coordinates().data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
  buffers_.push_back( points_buffer );

  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    gl_index& triangle_buffer = vao.triangles(k).buffer();

    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &triangle_buffer ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * vao.triangles(k).indices().size() , vao.triangles(k).indices().data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
    buffers_.push_back(triangle_buffer);
  }

  for (index_t k = 0; k < vao.nb_edges(); k++) {

    gl_index& edge_buffer = vao.edges(k).buffer();

    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &edge_buffer ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * vao.edges(k).indices().size() , vao.edges(k).indices().data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
    buffers_.push_back(edge_buffer);
  }

  for (index_t k = 0; k < vao.nb_fields(); k++) {

    // initialize the buffer
    gl_index& field_buffer = vao.field(k).buffer();
    gl_index& field_texture = vao.field(k).texture();

    // initialize the buffer
    GL_CALL( glGenBuffers( 1 , &field_buffer ) );
    GL_CALL( glGenTextures( 1 , &field_texture) );

    vao.field(k).write();

    buffers_.push_back(field_buffer);
    textures_.push_back(field_texture);
  }
}

OpenGL4_Manager::~OpenGL4_Manager() {
  glDeleteBuffers( buffers_.size() , buffers_.data() );
  glDeleteTextures( textures_.size() , textures_.data() );
  glDeleteVertexArrays( vao_.size() , vao_.data() );
}

void
PointPrimitive::draw( ShaderProgram& program ) {
  program.use();
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , buffer_) );
  GL_CALL( glDrawArrays( GL_POINTS , 0 , coordinates_.size()/3 ) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
}

void
EdgePrimitive::draw(bool with_tess) {

  // bind the buffer for the indices we want to draw
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );

  if (with_tess) {
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }
  else {
    GL_CALL( glDrawElements(GL_LINES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }
}

void
TrianglePrimitive::draw(bool with_tess) {
  // bind the buffer for the indices we want to draw
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
  if (with_tess) {
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }
  else {
    GL_CALL( glDrawElements(GL_TRIANGLES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }
}

void
FieldPrimitive::write() {

  FieldData* field = data_[active_name()].get();

  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , buffer_ ) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(gl_float) * field->data().size() , field->data().data() , GL_STATIC_DRAW) );

  GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) ); // fields are always in texture 0
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture_) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , buffer_ ) );
}

void
FieldPrimitive::activate( ShaderProgram& shader ) {

  // bind the desired solution texture to texture unit 0
  // solution fields are alway in texture unit 0 in the shaders
  glActiveTexture(GL_TEXTURE0 + 0);

  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture_) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , buffer_ ) );

  GLint solution_location = glGetUniformLocation(shader.handle() , "solution");
  glUniform1i(solution_location, 0); // first sampler in fragment shader
}

void
VertexAttributeObject::draw_triangles( ShaderProgram& shader ) {

  shader.use();
  shader.setUniform( "have_tessellation_shader" , shader.has_tessellation_shader() );

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  show_field_ = true;//false;
  if (solution_.size() == 0) show_field_ = false;

  // bind the desired colormap
  glActiveTexture(GL_TEXTURE0 + 1);
  GLint colormap_location = glGetUniformLocation(shader.handle() , "colormap");
  glUniform1i(colormap_location, 1); // second sampler in fragment shader

  for (index_t k = 0; k < triangles_.size(); k++) {
    if (!triangles_[k]->visible()) continue;

    if (show_field_) {
      solution_[k]->activate(shader);
      shader.setUniform( "use_constant_color" , 0 );
    }
    else {
      shader.setUniform( "use_constant_color" , 1 );
      shader.setUniform( "constant_color" , triangles_[k]->color() );
    }
    triangles_[k]->draw( shader.has_tessellation_shader() );
  }
}

void
VertexAttributeObject::draw_edges( ShaderProgram& shader ) {

  shader.use();
  shader.setUniform( "have_tessellation_shader" , shader.has_tessellation_shader() );

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  for (index_t k = 0; k < edges_.size(); k++) {
    edges_[k]->draw(shader.has_tessellation_shader());
  }
}

void
VertexAttributeObject::draw_points( ShaderProgram& shader ) {

  shader.use();

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  points_->draw(shader);
}

}

}
