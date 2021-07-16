#include "graphics/new/application.h"
#include "graphics/new/window.h"

#include "graphics/colormap.h"

namespace avro
{

namespace graphics
{

OpenGL_Application::OpenGL_Application() :
  window_(1024,768)
{
  window_.init();
  gui_ = std::make_shared<GUI>(window_);
}

void
OpenGL_Application::add( const TopologyBase& topology ) {
  std::shared_ptr<Plot> plot = std::make_shared<Plot>(topology);
  plot->build();
  plot_.push_back(plot);
  window_.add_plot(plot.get());
}

void
OpenGL_Application::run() {

  // bind the colormap values to a buffer
  gl_index colormap_buffer;
  GL_CALL( glGenBuffers( 1 , &colormap_buffer ) );
  Colormap colormap;
  colormap.change_style("giraffe");
  index_t ncolor = 256*3;
  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , colormap_buffer) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * ncolor , colormap.data() , GL_STATIC_DRAW) );

  // generate a texture to hold the colormap buffer
  GLuint colormap_texture;
  GL_CALL( glGenTextures( 1 , &colormap_texture) );
  GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , colormap_texture) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , colormap_buffer ) );

  // initial draw, subsequent drawing will only be performed when a callback is invoked
  window_.compute_view();
  gui_->draw();
  while (true) {

    // our thread will be put to sleep until user interaction is detected
    // so this does not actually redraw everything
    gui_->draw();

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window_.window())) break;
    if (glfwGetKey(window_.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}

} // graphics

} // avro
