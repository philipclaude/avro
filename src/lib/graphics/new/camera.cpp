#include "graphics/new/camera.h"

namespace avro
{

namespace graphics
{

Camera::Camera( float fov , index_t width , index_t height) :
  fov_(fov),
  width_(width),
  height_(height) {
  // compute the perspective matrix
  projection_matrix_ = glm::perspective(fov,float(width/height),0.01f,1000.0f);

  up_     = {0,1,0};
  lookat_ = {0,0,0};
  eye_    = {0,0,1};
}

void
Camera::compute_view() {
  view_matrix_ = glm::lookAt( eye_ , lookat_ , up_ );
}

}

}
