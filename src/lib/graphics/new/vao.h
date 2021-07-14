#ifndef AVRO_LIB_GRAPHICS_VAO_H_
#define AVRO_LIB_GRAPHICS_VAO_H_

#include "avro_types.h"

#include "graphics/gl.h"

#include <vector>

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

struct MeshFacet;

class VertexAttributeObject {

public:
  VertexAttributeObject( coord_t number , coord_t order ) :
    number_(number),
    order_(order),
    show_field_(true)
  {}

  void build( const TopologyBase& topology );

  void draw_triangles( ShaderProgram& );
  void draw_edges( ShaderProgram& );
  void draw_points( ShaderProgram& );

  void set_rank( index_t rank );
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

  std::shared_ptr<PointPrimitive> points_;
  std::vector< std::shared_ptr<EdgePrimitive> > edges_;
  std::vector< std::shared_ptr<TrianglePrimitive> > triangles_;
  std::vector< std::shared_ptr<FieldPrimitive> > solution_;

  mat4 model_matrix_;
};

} // graphics

} // avro

#endif
