#include "graphics/plotter.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"

#include "mesh/topology.h"

#include "numerics/field.h"

namespace ursa
{

namespace graphics
{

Primitive::Primitive( const TopologyHolder& topology , const Fields* fields ) :
  number_(topology.number()),
  topology_(topology),
  fields_(fields),
  active_(-1),
  shader_(NULL)
{}

void
Primitive::selectShader( Plotter* plotter )
{
  if (number_==2)
  {
    shader_ = &plotter->shader("basic");
    printf("selecting basic shader!\n");
  }
  else
    ursa_implement;

  // TODO: assign the shaders for the children too
}

ShaderProgram&
Primitive::shader()
{
  return *shader_;
}

void
WebGLPrimitive::write()
{
  // print to the web server

  // get the points

  // get the edges

  // get the triangles

  // send to wv
  handle_ = -1;

  ursa_implement;
}

void
WebGLPrimitive::draw()
{
  ursa_implement;
}

} // graphics

} // ursa
