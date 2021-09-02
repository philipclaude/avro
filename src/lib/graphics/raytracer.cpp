#include "graphics/gl.h"
#include "graphics/raytracer.h"
#include "graphics/shader.h"

#include "numerics/vec.hpp"

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
    data_[3*(i*width_+j)+d] = (*this)(i,j)[d];
}

void
Scene::intersect( const Ray& ray , Intersection& ixn ) const {



}

RayTracer::RayTracer( int width , int height ) :
  window_(width,height),
  canvas_(width,height)
{
  window_.init();
  canvas_.init_gl();

  window_.camera().set_lookat( {0.,0.,0.} );
  window_.camera().set_eye( {0.,0.,20.} );
  window_.camera().set_fov( M_PI/6.0 );

}

void
RayTracer::get_color( const Ray& ray , vec3& color ) {

  Intersection ixn;

  // determine which objects in the scene are hit
  scene_.intersect( ray , ixn );

  real_t t = 0.5*( ray.direction[1] ) + 0.2;

  vec3 cA = {1.,1.,1.};
  vec3 cB = {0.5,0.7,1.0};

  color = (1.0 - t)*cA + t*cB;
}

vec3
RayTracer::pixel2world( real_t u , real_t v ) const {

  // retrieve some parameters
  const real_t fov = window_.camera().fov();
  const real_t d = glm::norm( window_.camera().lookat() - window_.camera().eye() );
  const real_t a = real_t(canvas_.width()/canvas_.height());
  const real_t h = 2.0*d*tan(fov/2.0);
  const real_t w = a*h;

  // pixel coordiantes in camera space
  real_t pu = -0.5*w + w*u;
  real_t pv = -0.5*h + h*v;
  real_t pw = -d;
  vec3 q = {pu,pv,pw};

  // pixel coordinates in world space
  return basis_ * q + window_.camera().eye();
}

void
RayTracer::trace( index_t k ) {

  // get the pixel indices in [0,width] x [0,height]
  index_t i, j;
  pixel(k,i,j);

  // get the pixel coordinates in [0,1] x [0,1]
  real_t u = (j + 0.5)/window_.width();
  real_t v = (i + 0.5)/window_.height();

  // compute the world coordinates of the pixel
  vec3 p = pixel2world(u,v);

  // create a ray passing through the pixel
  Ray ray;
  ray.origin    = window_.camera().eye();
  ray.direction = glm::normalize( p - window_.camera().eye() );

  // determine which objects are intersected by the ray
  vec3 color;
  get_color( ray , color );

  canvas_(i,j) = color;
}

void
RayTracer::render() {

  const vec3& center = window_.camera().lookat();
  const vec3& eye = window_.camera().eye();
  vec3 up = {0.,1.,0.};

  // calculate the gaze
  vec3 g = center - eye;
  vec3 w = g * real_t(-1./glm::norm(g));
  vec3 u = glm::cross(up,w);
  u = u * (1./glm::norm(u));
  vec3 v = glm::cross(w,u);

  for (coord_t d = 0; d < 3; d++) {
    basis_(d,0) = u[d];
    basis_(d,1) = v[d];
    basis_(d,2) = w[d];
  }
  basis_.print();

  for (index_t k = 0; k < nb_pixels(); k++) {
    trace(k);
  }

  // determine if we want to render to the OpenGL framebuffer, or to an image
  if (true) {
    canvas_.convert();
    canvas_.draw_gl();
    glfwSwapBuffers(window_.window());
  }
}

} // graphics

} // avro
