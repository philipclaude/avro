#ifndef AVRO_LIB_GRAPHICS_VAO_H_
#define AVRO_LIB_GRAPHICS_VAO_H_

#include "avro_types.h"

#include "graphics/gl.h"

#include <vector>

#include <json/json.hpp>

namespace avro
{

class TopologyBase;
template<typename type> class Topology;

namespace graphics
{

class PointPrimitive;
class EdgePrimitive;
class TrianglePrimitive;
class FieldPrimitive;
class ClipPlane;

typedef struct
{
  coord_t dim;
  std::vector<index_t> indices;
} Facet;

struct CanonicalFacet : Facet
{
  index_t local;
};

struct MeshFacet : Facet
{
  std::vector<index_t> parent;
  std::vector<index_t> local;
  std::vector<int>     orientation;
};

class VertexAttributeObject {

public:
  VertexAttributeObject() :
    show_field_(true),
    uniform_color_(false),
    geometry_color_(false),
    tessellation_level_(8)
  {}

  void build( const TopologyBase& topology );

	void draw( const mat4& model , const mat4& view , const mat4& projection , const ClipPlane* clip=nullptr );
  void draw_triangles( ShaderProgram& );
  void draw_edges( ShaderProgram& );
  void draw_points( ShaderProgram& );

  void set_rank( index_t rank );
  void set_field( const std::string& name );
  bool& show_field() { return show_field_; }

  coord_t order() const { return order_; }

  index_t nb_edges() const { return edges_.size(); }
  index_t nb_triangles() const { return triangles_.size(); }
  index_t nb_fields() const { return solution_.size(); }

  const PointPrimitive& points() const { return *points_.get(); }
  const TrianglePrimitive& triangles( const index_t k ) const { return *triangles_[k].get(); }
  const EdgePrimitive& edges( index_t k ) const { return *edges_[k].get(); }
  const FieldPrimitive& field( index_t k ) const { return *solution_[k].get(); }

  PointPrimitive& points() { return *points_.get(); }
  TrianglePrimitive& triangles( const index_t k ) { return *triangles_[k].get(); }
  EdgePrimitive& edges( index_t k ) { return *edges_[k].get(); }
  FieldPrimitive& field( index_t k ) { return *solution_[k].get(); }

  gl_index& vertex_array()       { return vertex_array_; }
  gl_index  vertex_array() const { return vertex_array_; }

  void apply_transformation( const mat4& m );

  const nlohmann::json& get_info() const { return info_; }
  int& tessellation_level() { return tessellation_level_; }
  bool& uniform_color() { return uniform_color_; }
  bool& geometry_color() { return geometry_color_; }

private:
  template<typename type>
  void get_primitives( const Topology<type>& topology , const std::vector<std::vector<MeshFacet>>& facets );

  template<typename type>
  void _build( const Topology<type>& topology );

private:
  coord_t number_;
  coord_t order_;

  gl_index vertex_array_;
  bool show_field_;
  bool uniform_color_;
  bool geometry_color_;

  std::shared_ptr<PointPrimitive> points_;
  std::vector< std::shared_ptr<EdgePrimitive> > edges_;
  std::vector< std::shared_ptr<TrianglePrimitive> > triangles_;
  std::vector< std::shared_ptr<FieldPrimitive> > solution_;

  mat4 model_matrix_;

  nlohmann::json info_;
  int tessellation_level_;
};

} // graphics

} // avro

#endif
