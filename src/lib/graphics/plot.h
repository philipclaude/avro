#ifndef avro_LIB_GRAPHICS_PLOT_H_
#define avro_LIB_GRAPHICS_PLOT_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <memory>
#include <vector>

namespace avro
{

class Fields;

namespace graphics
{

class Controls;
class Plotter;
class Primitive;
class Window;

class Plot
{
public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  Plot( const TopologyBase& topology , Window* window );

  void draw();
  void write();

  void add( Primitive_ptr prim );
  void setWindow( Window* window );

  void set_transform_feedback( bool x );

  const TopologyBase& topology() const { return topology_; }

  Primitive& primitive( index_t k );
  void set_visibility( const Controls& controls );

private:
  const TopologyBase& topology_;

  std::vector< Primitive_ptr > primitive_;
  Window* window_;
  Plotter* plotter_;
};



} // graphics

} // avro

#endif
