#include "graphics/colormap.h"
#include "graphics/application.h"
#include "graphics/managers.h"
#include "graphics/plot.h"
#include "graphics/primitives.h"
#include "graphics/shader_library.h"
#include "graphics/vao.h"

namespace avro
{

namespace graphics
{

std::shared_ptr<Shaders> __shaders__;

OpenGL4_Manager::OpenGL4_Manager() {

  // initialize GLFW
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  glfwWindowHint(GLFW_RESIZABLE , GLFW_TRUE );
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
  vao.points().memory() = sizeof(gl_float) * vao.points().coordinates().size();

  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    gl_index& triangle_buffer = vao.triangles(k).buffer();

    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &triangle_buffer ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_buffer ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * vao.triangles(k).indices().size() , vao.triangles(k).indices().data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
    buffers_.push_back(triangle_buffer);
    vao.triangles(k).memory() = sizeof(gl_index) * vao.triangles(k).indices().size();
  }

  for (index_t k = 0; k < vao.nb_edges(); k++) {

    gl_index& edge_buffer = vao.edges(k).buffer();

    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &edge_buffer ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * vao.edges(k).indices().size() , vao.edges(k).indices().data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
    buffers_.push_back(edge_buffer);
    vao.edges(k).memory() = sizeof(gl_index) * vao.edges(k).indices().size();
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
  glfwTerminate();
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

  memory_ = sizeof(gl_float) * field->data().size();
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
  shader.setUniform( "u_level" , tessellation_level_ );

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  show_field_ = true;
  if (solution_.size() == 0) show_field_ = false;

  // bind the desired colormap
  glActiveTexture(GL_TEXTURE0 + 1);
  GLint colormap_location = glGetUniformLocation(shader.handle() , "colormap");
  glUniform1i(colormap_location, 1); // second sampler in fragment shader

  for (index_t k = 0; k < triangles_.size(); k++) {
    if (!triangles_[k]->visible()) continue;

    shader.setUniform("u_alpha",float(1.0));
    if (show_field_) {
      solution_[k]->activate(shader);
      shader.setUniform( "use_constant_color" , 0 );

      const std::vector<gl_float>& data = solution_[k]->active().data();
      float umin = * std::min_element( data.begin() , data.end() );
      float umax = * std::max_element( data.begin() , data.end() );

      shader.setUniform( "u_umin" , umin );
      shader.setUniform( "u_umax" , umax );
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
  shader.setUniform( "u_level" , tessellation_level_ );

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

void
VertexAttributeObject::draw( const mat4& model , const mat4& view , const mat4& projection , const ClipPlane* clip ) {

  if (number_ == 2) {
    glDisable(GL_CULL_FACE);
  }
  else {
    glEnable(GL_CULL_FACE);
  }

  // bind which attributes we want to draw
  GL_CALL( glBindVertexArray(vertex_array_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_->buffer()  ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

  // calculate the matrices for the shaders
  mat4 mv = view * model;
  mat4 mvp = projection * mv;
  mat4 normal_matrix = glm::transpose(glm::inverse(mv));

  // draw the edges
  for (index_t k = 0; k < edges_.size(); k++) {

    if (!edges_[k]->visible()) continue;
    coord_t q = edges_[k]->order();

    ShaderProgram& shader = __shaders__->get("edges",0,q, q > 1);
    shader.use();
    shader.setUniform("u_ModelViewProjectionMatrix",mvp);
    if (shader.has_tessellation_shader())
      shader.setUniform( "u_level" , tessellation_level_ );

    if (clip != nullptr && clip->style() > 0) {
      vec3 normal;
      vec3 center;
      clip->get(center,normal);

      shader.setUniform("u_clip",clip->style());
      shader.setUniform("u_clip_center",center);
      shader.setUniform("u_clip_normal",normal);
    }
    else {
      shader.setUniform("u_clip",-1);
    }

    edges_[k]->draw(shader.has_tessellation_shader());
  }

  // draw the triangles
  show_field_ = (!uniform_color_ && !geometry_color_);
  if (solution_.size() == 0) show_field_ = false;
  for (index_t k = 0; k < triangles_.size(); k++) {

    if (!triangles_[k]->visible()) continue;

    coord_t p = 0;
    if (solution_.size() != 0) p = solution_[k]->active().order();
    coord_t q = triangles_[k]->order();

    ShaderProgram& shader = (show_field_) ? __shaders__->get("triangles",p,q, q > 1) : __shaders__->get("triangles",-1,q,q>1);
    //ShaderProgram& shader = __shaders__->get("triangles",p,q, q > 1);
    shader.use();
    shader.setUniform("u_ModelViewProjectionMatrix",mvp);
    shader.setUniform("u_NormalMatrix",normal_matrix);
    shader.setUniform("u_ModelViewMatrix",mv);
    shader.setUniform("u_alpha",float(1.0));

    if (clip != nullptr && clip->style() > 0) {
      vec3 normal;
      vec3 center;
      clip->get(center,normal);

      shader.setUniform("u_clip",clip->style());
      shader.setUniform("u_clip_center",center);
      shader.setUniform("u_clip_normal",normal);
    }
    else {
      shader.setUniform("u_clip",-1);
    }

    if (shader.has_tessellation_shader()) {
      shader.setUniform( "u_level" , tessellation_level_ );
    }

    // bind the desired colormap
    glActiveTexture(GL_TEXTURE0 + 1);
    GLint colormap_location = glGetUniformLocation(shader.handle() , "colormap");
    glUniform1i(colormap_location, 1); // second sampler in fragment shader

    if (number_ == 2) shader.setUniform("u_lighting",-1);
    else {
      if (lighting_)
        shader.setUniform("u_lighting",1);
      else
        shader.setUniform("u_lighting",-1);
    }

    //shader.setUniform("u_lighting",-1);

    if (show_field_) {
      solution_[k]->activate(shader);
      shader.setUniform( "use_constant_color" , 0 );

      const std::vector<gl_float>& data = solution_[k]->active().data();
      float umin = * std::min_element( data.begin() , data.end() );
      float umax = * std::max_element( data.begin() , data.end() );

      shader.setUniform( "u_umin" , umin );
      shader.setUniform( "u_umax" , umax );

    }
    else if (geometry_color_) {
      shader.setUniform( "use_constant_color" , 1 );
      shader.setUniform( "constant_color" , triangles_[k]->color() );
    }
    else {
      vec3 c = {1.0,1.0,0.8}; // yellow
      shader.setUniform( "use_constant_color" , 1 );
      shader.setUniform( "constant_color" , c );
    }
    triangles_[k]->draw( shader.has_tessellation_shader() );
  }

  // draw the points
  // TODO
}

void
ClipPlane::write() {

  GL_CALL( glGenVertexArrays( 1, &vertex_array_ ) );
  GL_CALL( glBindVertexArray(vertex_array_) );

  gl_index buffers[2];
  GL_CALL( glGenBuffers( 2 , buffers ) );
  point_buffer_ = buffers[0];
  index_buffer_ = buffers[1];

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, point_buffer_ ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * 12 , coordinates_ , GL_STATIC_DRAW) );

  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_ ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * 6 , indices_ , GL_STATIC_DRAW) );
}

void
ClipPlane::draw( const mat4& view , const mat4& projection ) const {

  if (!visible_) return;

  // calculate the matrices for the shaders
  mat4 mv = view * plot_.model_matrix() * transform_matrix_;
  mat4 mvp = projection * mv;

  // pick an appropriate program
  ShaderProgram& program = __shaders__->get("triangles",0,1,false); // p = 0, q = 1, no tessellation shader
  program.use();

  vec3 color = {0.9,0.2,0.2}; // pink
  program.setUniform("use_constant_color",1);
  program.setUniform("constant_color",color);
  program.setUniform("u_alpha",float(0.2));
  program.setUniform("u_lighting",0); // no lighting
  program.setUniform("u_clip",-1);

  program.setUniform("u_ModelViewProjectionMatrix",mvp);
  program.setUniform("u_ModelViewMatrix",mv);

  GL_CALL( glBindVertexArray(vertex_array_) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, point_buffer_ ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );

  glDisable(GL_CULL_FACE);
  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_buffer_ ) );
  GL_CALL( glDrawElements(GL_TRIANGLES, 6 , GL_UNSIGNED_INT , 0 ) );
}

ClipPlane::~ClipPlane() {
  glDeleteBuffers(1,&point_buffer_);
  glDeleteBuffers(1,&index_buffer_);
  glDeleteVertexArrays(1,&vertex_array_);
}

void
Window::select_colormap(const std::string& name) {

  // bind the colormap values to a buffer
  gl_index colormap_buffer;
  GL_CALL( glGenBuffers( 1 , &colormap_buffer ) );
  Colormap colormap;
  colormap.change_style(name);
  index_t ncolor = 256*3;
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , colormap_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * ncolor , colormap.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the colormap buffer
  GLuint colormap_texture;
  GL_CALL( glGenTextures( 1 , &colormap_texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , colormap_texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , colormap_buffer ) );

  manager_.track_texture(colormap_texture);
  manager_.track_buffer(colormap_buffer);
}

}

}
