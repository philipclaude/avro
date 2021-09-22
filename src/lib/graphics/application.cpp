//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/parallel_for.h"
#include "graphics/application.h"
#include "graphics/camera.h"
#include "graphics/window.h"

#include "graphics/colormap.h"

#include <json/json.hpp>

#include <fstream>
#include <unistd.h>

namespace avro
{

namespace graphics
{

OpenGL_Application::OpenGL_Application() :
  //window_(1024,768)
  window_(800,800)
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
OpenGL_Application::run(bool quit) {

  #if AVRO_WITH_GL

  if (quit)
    return;

  // initial draw, subsequent drawing will only be performed when a callback is invoked
  window_.compute_view();
  gui_->draw();

  while (true) {

    // our thread will be put to sleep until user interaction is detected
    // so this does not actually redraw everything, unless interaction is detected
    gui_->draw();

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window_.window())) break;
    if (glfwGetKey(window_.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
  #else
  avro_assert_not_reached;
  #endif
}

void
WebGL_Application::add( const TopologyBase& topology ) {
  std::shared_ptr<Plot> plot = std::make_shared<Plot>(topology,false);
  plot->build();
  plot_.push_back(plot);
}

class RunThread
{
public:
  typedef RunThread thisclass;

  RunThread(WebGL_Application& app,bool quit) :
    app_(app),
    quit_(quit) {}

  void run( index_t i )
  {
    if (i == 0) {
      app_.run_thread(quit_);
    }
    else if (i == 1) {
      //usleep(1000);
      //std::string cmd = "open " + AVRO_SOURCE_DIR + "/app/avro.html";
      //system(cmd.c_str());
    }
  }

  void runapp( const index_t nthread )
  {
    ProcessCPU::parallel_for (
      parallel_for_member_callback( this , &thisclass::run ),
      0,nthread
    );
  }
private:
  WebGL_Application& app_;
  bool quit_;
};

void
WebGL_Application::run_thread(bool quit) {
  if (quit) return;
  manager_.send(7681);
}

void
WebGL_Application::run(bool quit) {

  for (index_t k = 0; k < plot_.size(); k++) {
    for (index_t j = 0; j < plot_[k]->nb_vao(); j++)
      manager_.write( plot_[k]->vao(j) );
  }

  // we will start up two threads, one to write, and one to call the browser (after waiting a bit)
  RunThread runner(*this,quit);
  runner.runapp(2);
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
Viewer::run(bool quit) {
  app_->run(quit);
}

} // graphics

} // avro
