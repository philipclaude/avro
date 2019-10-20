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
  printf("topology number = %u\n",topology_.number());
  Primitive_ptr prim = std::make_shared<OpenGLPrimitive>(topology,window);
  add(prim);

  // assign the shader for this primitive
  prim->selectShader( plotter_ );

  prim->shader().use();
  prim->shader().printActiveAttribs();
  prim->shader().printActiveUniforms();

  prim->shader().setUniform( "MVP" , window->mvp() );
  //prim->shader().setUniform( "u_modelViewMatrix" , window->mv() );
  prim->shader().setUniform( "u_normalMatrix" , window->normal() );
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


} // graphics

} // ursa
