#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include "numerics/field.h"

namespace ursa
{

namespace graphics
{

Plot::Plot( const TopologyHolder& topology , Fields* fields ) :
  topology_(topology),
  fields_(fields),
  window_(NULL),
  plotter_(NULL)
{
  printf("topology number = %u\n",topology_.number());
  Primitive_ptr prim = std::make_shared<OpenGLPrimitive>(topology);
  add(prim);

  // assign the shader for this primitive
  prim->selectShader( plotter_ );
}

void
Plot::setWindow( Window* window )
{
  window_  = window;
  plotter_ = window->plotter();
}

void
Plot::draw()
{
  printf("drawing!\n");

  // loop through the primitives
  for (index_t k=0;k<primitive_.size();k++)
  {
    //primitive_[k]->shader().setUniforms(*window_);
    primitive_[k]->draw();
  }
}

void
Plot::write()
{
  printf("writing!\n");
  // loop through the primitives
  for (index_t k=0;k<primitive_.size();k++)
    primitive_[k]->write();
}

void
Plot::add( Primitive_ptr prim )
{
  primitive_.push_back( prim );
}

DummyTopology::DummyTopology() :
  Topology<Simplex<Lagrange>>(vertices_,2),
  vertices_(3)
{
  real_t x0[3] = {0,0,0};
  real_t x1[3] = {1,0,0};
  real_t x2[3] = {0,1,0};
  real_t x3[3] = {1,1,0};
  vertices_.create( x0 );
  vertices_.create( x1 );
  vertices_.create( x2 );
  vertices_.create( x3 );

  index_t t0[3] = {0,1,2};
  index_t t1[3] = {0,3,2};

  add( t0 , 3 );
  add( t1 , 3 );
}



} // graphics

} // ursa
