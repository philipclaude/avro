#include "graphics/gl.h"
#include "graphics/raytracer.h"
#include "graphics/shader.h"

namespace avro
{

namespace graphics
{

Canvas::Canvas( int width , int height ) :
  width_(width),
  height_(height),
  pixel_(width_*height_)
{}

void
Canvas::init_gl() {

  // only do the following if opengl is supported
  GL_CALL( glGenVertexArrays( 1, &vertex_array_ ) );
  GL_CALL( glBindVertexArray(vertex_array_) );

  // bind the colormap values to a buffer
  GL_CALL( glGenBuffers( 1 , &pixel_buffer_ ) );
  GL_CALL( glGenTextures( 1 , &pixel_texture_ ) );

  // initialize the shader
  std::vector<std::string> macros = {"#version 410"};
  shader_ = std::make_shared<ShaderProgram>("raytracer",false,macros);
  shader_->use();

  shader_->setUniform("u_width",width_);
  shader_->setUniform("u_height",height_);

  // bind the desired texture
  glActiveTexture(GL_TEXTURE0 + 0);
  GLint pixels_location = glGetUniformLocation(shader_->handle() , "pixels");
  glUniform1i(pixels_location, 0); // first sampler in fragment shader
}

void
Canvas::draw_gl() {

  // clear the screen
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  // write the pixel data to the pixel buffer
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , pixel_buffer_) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(gl_float)*data_.size() , data_.data() , GL_STATIC_DRAW) );

  // bind the pixel buffer to the pixel texture
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , pixel_texture_ ) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_RGB32F , pixel_buffer_ ) );

  // nothing is actually drawn here, we just rely on the interpolation to give (u,v) coordinates to then look up the texture value
  GL_CALL( glDrawArrays( GL_POINTS , 0 , 1 ) );
}

void
Canvas::convert() {

  data_.resize( width_*height_*3 , 0.0 );

  for (index_t i = 0; i < height_; i++)
  for (index_t j = 0; j < width_; j++)
  for (index_t d = 0; d < 3; d++)
    data_[3*(i*width_+j)+d] = real_t(j)/width_;

  print_inline(data_);
}

RayTracer::RayTracer( int width , int height ) :
  window_(width,height),
  canvas_(width,height)
{
  window_.init();
  canvas_.init_gl();
}

void
RayTracer::draw() {

  printf("draw scene!\n");

  // determine if we want to render to the OpenGL framebuffer, or to an image
  if (true) {

    canvas_.convert();
    canvas_.draw_gl();

    glfwSwapBuffers(window_.window());

  }
}



} // graphics

} // avro
