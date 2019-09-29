#include "graphics/primitive.h"

#include "numerics/field.h"

namespace ursa
{

namespace graphics
{

Primitive::Primitive( const TopologyHolder& topology , const Fields* fields ) :
  topology_(topology),
  fields_(fields),
  identifier_(-1),
  active_(-1)
{}

void
WebGLPrimitive::write()
{
  // print to the web server

  // get the points

  // get the edges

  // get the triangles

  // send to wv
  identifier_ = -1;

}

void
OpenGLPrimitive::write()
{
  // bind the buffers to the opengl context


}

} // graphics

} // ursa
