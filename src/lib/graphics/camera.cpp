#include "graphics/camera.h"

namespace avro
{

namespace graphics
{

Camera::Camera( float fov , index_t width , index_t height) :
  fov_(fov),
  width_(width),
  height_(height) {
  up_     = {0,1,0};
  lookat_ = {0.5,0.5,0.5};
  eye_    = {0,0,5};
  compute_projection(width,height);
  compute_view();
}

void
Camera::compute_view() {
  view_matrix_ = glm::lookAt( eye_ , lookat_ , up_ );
}

void
Camera::compute_projection( index_t width , index_t height ) {
  width_  = width;
  height_ = height;
  projection_matrix_ = glm::perspective(fov_,float(width)/float(height),0.001f,1000.0f);
}

}

}
