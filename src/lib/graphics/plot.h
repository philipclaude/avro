#ifndef URSA_LIB_GRAPHICS_PLOT_H_
#define URSA_LIB_GRAPHICS_PLOT_H_

#include <memory>
#include <vector>

namespace ursa
{

namespace graphics
{

class Primitive;
class Fields;
class TopologyHolder;

class Plot
{
public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  Plot( const TopologyHolder& topology , const Fields* fields );

private:
  const TopologyHolder& topology_;
  const Fields* fields_;

  std::vector< Primitive_ptr > primitive_;
};

} // graphics

} // ursa

#endif
