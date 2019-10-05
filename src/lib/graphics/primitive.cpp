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
    shader_ = &plotter->shader("wv");
    printf("selecting wv shader!\n");

    // set some uniforms
    shader_->use();
    shader_->setUniform( "lightDir" , 0.0 , 0.3 , 1.0 );
    shader_->setUniform( "wAmbient" , 0.25f );
    shader_->setUniform( "xpar" , 1.0f );
    shader_->setUniform( "conNormal" , 0. , 0. , 1. );
    shader_->setUniform( "conColor" , 1. , 1. , 0. );
    shader_->setUniform( "bacColor" , 0.5 , 0.5 , 0.5 );
    shader_->setUniform( "wColor" , 0.0f );
    shader_->setUniform( "bColor" , 1.0f );
    shader_->setUniform( "wNormal" , 0.0f );
    shader_->setUniform( "wLight" , 1.0f );
    shader_->setUniform( "pointSize" , 2.0f );
    shader_->setUniform( "picking" , 0 );
    shader_->setUniform( "vbonum" , 0 );
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
