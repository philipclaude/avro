#ifndef AVRO_LIB_GRAPHICS_MESH_H_
#define AVRO_LIB_GRAPHICS_MESH_H_

#include "avro_types.h"

#include <vector>

namespace avro
{

// forward declarations
class TopologyBase;
template<typename type> class Topology;

namespace graphics
{
struct MeshFacet;

class Tessellation {

public:
  Tessellation( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  void build( const TopologyBase& topology );

private:
  void get_primitives( const TopologyBase& topology , const std::vector<std::vector<MeshFacet>>& facets );

  template<typename type>
  void _build( const Topology<type>& topology );

private:
  coord_t number_;
  coord_t order_;

  /*PointPrimitive points_;
  std::vector<EdgePrimitive> edges_;
  std::vector<TrianglePrimitive> triangles_;
  std::vector<SolutionPrimitive> solution_;*/
};

} // graphics

} // avro

#endif
