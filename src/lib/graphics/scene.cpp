#include "graphics/controls.h"
#include "graphics/scene.h"

namespace avro
{

namespace graphics
{

void
SceneGraph::update_matrices( const Controls& controls )
{
  normal_matrix_ = controls.normal();
  mvp_matrix_ = controls.model_view_projection();
}

} // graphics

} // avro
