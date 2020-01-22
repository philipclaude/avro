#ifndef avro_LIB_GRAPHICS_SCENE_H_
#define avro_LIB_GRAPHICS_SCENE_H_

#include "graphics/controls.h"
#include "graphics/math.h"
#include "graphics/primitive.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace avro
{

namespace graphics
{

class Primitive;

class SceneGraph
{
public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  SceneGraph() :
    update_(true)
  {}

  index_t add_primitive( const TopologyBase& topology )
  {
    index_t id = primitive_.size();
    Primitive_ptr primitive = std::make_shared<Primitive>(topology,this);
    primitive_.push_back(primitive);
    return id;
  }

  void write( GraphicsManager& manager )
  {
    for (index_t k=0;k<primitive_.size();k++)
      primitive_[k]->write(manager);
  }

  bool update() const { return update_; }
  void set_update( bool x ) { update_ = x; }

  void update_matrices( const Trackball& trackball , float,float,float );

  const mat4& mvp_matrix() const { return mvp_matrix_; }
  const mat4& normal_matrix() const { return normal_matrix_; }

  mat4& mvp_matrix() { return mvp_matrix_; }
  mat4& normal_matrix() { return normal_matrix_; }

  index_t nb_primitives() const { return primitive_.size(); }

  Primitive& primitive( index_t k ) { return *primitive_[k].get(); }

private:
  std::vector<Primitive_ptr> primitive_; // roots of the scene graph

  // store all the matrices here
  mat4 mvp_matrix_;
  mat4 view_matrix_;
  mat4 proj_matrix_;
  mat4 model_matrix_;
  mat4 normal_matrix_;

  bool update_;
};

} // graphics

} // avro

#endif
