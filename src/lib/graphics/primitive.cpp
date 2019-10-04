#include "graphics/plotter.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"

#include "mesh/field.h"
#include "mesh/topology.h"

namespace ursa
{

namespace graphics
{

Primitive::Primitive( const TopologyHolder& topology , Window* window ) :
  number_(topology.number()),
  topology_(topology),
  active_(-1),
  shader_(NULL),
  window_(window)
{}

void
Primitive::selectShader( Plotter* plotter )
{
  if (number_==2)
  {
    shader_ = &plotter->shader("basic");
    printf("selecting basic shader!\n");
  }
  else if (number_==1)
  {
    shader_ = &plotter->shader("edge");
    printf("selecting edge shader!\n");
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
