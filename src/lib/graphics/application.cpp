#include "graphics/application.h"

namespace avro
{

namespace graphics
{

void
ApplicationBase::write()
{
  // write all the data present in the scene graph
  for (index_t k=0;k<scene_.size();k++)
    scene_[k].write(manager_);
}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;
template class Application<Web_Interface>;

} // graphics

} // avro
