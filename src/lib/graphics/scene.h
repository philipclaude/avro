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

  SceneGraph() :
    update_(true)
  {
    menu_["primitives"] = {};
  }

  index_t add_primitive( const TopologyBase& topology )
  {
    index_t id = primitive_.size();
    Primitive_ptr primitive = std::make_shared<Primitive>(topology,this);
    primitive_.push_back(primitive);

    std::vector<const TopologyBase*> children;
    topology.get_topologies(children);
    for (index_t k=0;k<children.size();k++)
    {
      //if (children[k]->number()<1) continue;
      primitive->add_child( std::make_shared<Primitive>(*children[k],this) );
    }
    return id;
  }

  const json& menu() const { return menu_; }

  void write( GraphicsManager& manager )
  {
    for (index_t k=0;k<primitive_.size();k++)
      primitive_[k]->write(manager);

    std::vector<std::string> primitives;
    for (index_t k=0;k<primitive_.size();k++)
    {
      MenuEntry entry(*primitive_[k].get());
      primitives.push_back( entry.dump() );
    }
    menu_["primitives"] = primitives;
  }

  void remove( index_t k )
  {
    primitive_.erase( primitive_.begin()+k );
  }

  bool update() const { return update_; }
  void set_update( bool x ) { update_ = x; }

  void update_matrices( const Controls& controls );

  const mat4& mvp_matrix() const { return mvp_matrix_; }
  const mat4& normal_matrix() const { return normal_matrix_; }

  mat4& mvp_matrix() { return mvp_matrix_; }
  mat4& normal_matrix() { return normal_matrix_; }

  index_t nb_primitives() const { return primitive_.size(); }

  Primitive& primitive( index_t k ) { return *primitive_[k].get(); }

private:
  std::vector<Primitive_ptr> primitive_; // roots of the scene graph
  json menu_;

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
