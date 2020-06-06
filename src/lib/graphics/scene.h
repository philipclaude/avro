#ifndef avro_LIB_GRAPHICS_SCENE_H_
#define avro_LIB_GRAPHICS_SCENE_H_

#include "common/json.h"

#include "graphics/math.h"
#include "graphics/primitive.h"

#include "mesh/topology.h"

#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <json/json.hpp>

namespace avro
{

namespace graphics
{

class Primitive;
class Controls;

class SceneGraph
{
private:
  class MenuEntry : public json
  {
  public:
    MenuEntry( const Primitive& primitive )
    {
      const void* address = static_cast<const void*>(&primitive);
      std::stringstream ss;
      ss << address;
      this->operator[]("name") = ss.str()+"-root";
      add(primitive);

      this->operator[]("Volumes") = volumes_;
      this->operator[]("Faces") = faces_;
      this->operator[]("Edges") = edges_;
      this->operator[]("Nodes") = nodes_;

      std::vector<std::string> field_names;
      primitive.topology().fields().get_names(field_names);
      this->operator[]("fields") = field_names;
    }

    void add( const Primitive& primitive )
    {
      // check if this primitive has been added
      if (primitives_.find(&primitive)!=primitives_.end()) return;
      primitives_.insert( &primitive );

      const void* address = static_cast<const void*>(&primitive);
      std::stringstream ss;
      ss << address;

      if (primitive.number()==4)
      {
        avro_implement; // extract boundary
      }

      if (primitive.number()==3)
      {
        volumes_.push_back( ss.str() );
      }
      if (primitive.number()==2)
      {
        faces_.push_back( ss.str() );
      }
      if (primitive.number()==1)
      {
        edges_.push_back( ss.str() );
      }
      if (primitive.number()==0)
      {
        nodes_.push_back( ss.str() );
      }

      for (index_t k=0;k<primitive.nb_children();k++)
        add(primitive.child(k));
    }

  private:
    std::set<const Primitive*> primitives_;

    std::vector<std::string> volumes_;
    std::vector<std::string> faces_;
    std::vector<std::string> edges_;
    std::vector<std::string> nodes_;
  };

public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  SceneGraph();

  index_t add_primitive( const TopologyBase& topology );
  void remove( index_t k );

  const json& menu() const { return menu_; }

  void write( GraphicsManager& manager );

  bool update() const { return update_; }
  void set_update( bool x ) { update_ = x; }

  void update_matrices( const Controls& controls );

  const mat4& mvp_matrix() const { return mvp_matrix_; }
  const mat4& normal_matrix() const { return normal_matrix_; }

  mat4& mvp_matrix() { return mvp_matrix_; }
  mat4& normal_matrix() { return normal_matrix_; }

  index_t nb_primitives() const { return primitive_.size(); }

  Primitive& primitive( index_t k ) { return *primitive_[k].get(); }

  void get_bounding_box( real_t* box ) const;
  void set_focus( real_t* focus );

  const real_t* focus() const { return focus_; }

private:
  std::vector<Primitive_ptr> primitive_; // roots of the scene graph
  json menu_;

  // store all the matrices here
  mat4 mvp_matrix_;
  mat4 normal_matrix_;

  bool update_;

  real_t focus_[4];
};

} // graphics

} // avro

#endif
