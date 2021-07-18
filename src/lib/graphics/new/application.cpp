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

void
WebGL_Application::add( const TopologyBase& topology ) {
  std::shared_ptr<Plot> plot = std::make_shared<Plot>(topology,false);
  plot->build();
  plot_.push_back(plot);
}

void
WebGL_Application::run() {

  for (index_t k = 0; k < plot_.size(); k++) {
    for (index_t j = 0; j < plot_[k]->nb_vao(); j++)
      manager_.write( plot_[k]->vao(j) );
  }

  manager_.send(7681);
}

Viewer::Viewer(bool web) {

  if (!web) {
    try {
      app_ = std::make_shared<OpenGL_Application>();
    }
    catch (...) {
      printf("[warning] it doesn't seem OpenGL4 is supported, using the web viewer");
      web = true;
    }
  }

  if (web) app_ = std::make_shared<WebGL_Application>();
}

void
Viewer::add( const TopologyBase& topology ) {
  app_->add(topology);
}

void
Viewer::run() {
  app_->run();
}

} // graphics

} // avro
