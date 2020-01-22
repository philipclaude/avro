#include "graphics/scene.h"

namespace avro
{

namespace graphics
{

void
SceneGraph::update_matrices( const Trackball& trackball , float fov , float width , float height )
{
  model_matrix_ = mat4(1.0);
  proj_matrix_  = glm::perspective(glm::radians(fov), float(width)/float(height) , 0.1f, 100.0f);

  // compute the matrices that need to be passed to the shaders
  view_matrix_   = trackball.camera().view_matrix;
  mvp_matrix_    = proj_matrix_ * view_matrix_ * model_matrix_;
  normal_matrix_ = glm::transpose(glm::inverse(glm::mat3( model_matrix_*view_matrix_)));
}

} // graphics

} // avro
