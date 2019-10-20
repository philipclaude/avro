#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include "mesh/field.h"

namespace ursa
{

namespace graphics
{

Plot::Plot( const TopologyHolder& topology , Window* window ) :
  topology_(topology),
  window_(window),
  plotter_(window->plotter())
{
  Primitive_ptr prim = std::make_shared<OpenGLPrimitive>(topology,window);
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
  // loop through the primitives
  for (index_t k=0;k<primitive_.size();k++)
  {
    primitive_[k]->draw();
  }
}

void
Plot::write()
{
  // loop through the primitives
  for (index_t k=0;k<primitive_.size();k++)
    primitive_[k]->write();
}

Primitive&
Plot::primitive( index_t k )
{
  return *(primitive_[k].get());
}

void
Plot::add( Primitive_ptr prim )
{
  primitive_.push_back( prim );
}

} // graphics

} // ursa
