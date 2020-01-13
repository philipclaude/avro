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

void
OpenGL_Manager::draw()
{
  // draws with the current shader (selected by the scene graph)
  std::map<index_t,index_t>::iterator it;

  for (it=vao_triangles_.begin();it!=vao_triangles_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_triangles_active_[vao];
    ShaderProgram* shader = vao_triangles_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glDrawElements(GL_TRIANGLES,nb, GL_UNSIGNED_INT , 0 ) )
  }

  for (it=vao_edges_.begin();it!=vao_edges_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_edges_active_[vao];
    ShaderProgram* shader = vao_edges_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glDrawElements( GL_LINES , nb , GL_UNSIGNED_INT , 0 ) )
  }

  for (it=vao_points_.begin();it!=vao_points_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_points_active_[vao];
    ShaderProgram* shader = vao_points_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , nb ) );
  }
}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

} // graphics

} // avro
