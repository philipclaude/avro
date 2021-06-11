//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/application.h"
#include "graphics/gl.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include <wsserver.h>

#include <limits>

namespace avro
{

namespace graphics
{

void
ApplicationBase::write()
{
  // write all the data present in the scene graph
  for (index_t k=0;k<scenes_.size();k++)
    scenes_[k]->write(manager_);
}

void
ApplicationBase::focus_scenes()
{
  // initialize the bounding box to something really far away
  bounding_box_[0] = bounding_box_[1] = bounding_box_[2] = std::numeric_limits<real_t>::max();
  bounding_box_[3] = bounding_box_[4] = bounding_box_[5] = std::numeric_limits<real_t>::min();

  for (index_t k=0;k<scenes_.size();k++)
  {
    // update the bounding box based on the points in each scene topology
    scenes_[k]->get_bounding_box(bounding_box_);
  }

  printf("--> bounding box: x = (%g,%g), y = (%g,%g), z = (%g,%g)\n",
              bounding_box_[0],bounding_box_[3],
              bounding_box_[1],bounding_box_[4],
              bounding_box_[2],bounding_box_[5]);


  // calculate the focus
  focus_[0] = 0.5*( bounding_box_[0] + bounding_box_[3] );
  focus_[1] = 0.5*( bounding_box_[1] + bounding_box_[4] );
  focus_[2] = 0.5*( bounding_box_[2] + bounding_box_[5] );
  focus_[3] = std::sqrt( (bounding_box_[3] - bounding_box_[0])*(bounding_box_[3] - bounding_box_[0])
                        +(bounding_box_[4] - bounding_box_[1])*(bounding_box_[4] - bounding_box_[1])
                        +(bounding_box_[5] - bounding_box_[2])*(bounding_box_[5] - bounding_box_[2]) );

  printf("--> focus: center = (%g,%g,%g), scale = %g\n",focus_[0],focus_[1],focus_[2],focus_[3]);

  // focus all the scenes
  for (index_t k=0;k<scenes_.size();k++)
  {
    // update the bounding box based on the points in each scene topology
    scenes_[k]->set_focus(focus_);
  }
}

void
Application<Web_Interface>::send( const std::string& response ) const
{
  wv_broadcastText( const_cast<char*>(response.c_str()) );
}

void
Application<Web_Interface>::save_eps()
{
  #ifdef AVRO_WITH_GL
  // first set the transformation to the scene
  OpenGL_Manager manager_gl;
  scene_->write(manager_gl);
  TransformFeedbackResult feedback;
  manager_gl.draw(*scene_.get(),&feedback);
  #else
  avro_assert_not_reached;
  #endif
}

#ifdef AVRO_WITH_GL

template<typename API_t>
Application<GLFW_Interface<API_t>>::Application() :
  ApplicationBase(manager_),
  restart_(false)
{
  // initialize OpenGL
  avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );

  // set the version
  #ifdef AVRO_HEADLESS_GRAPHICS // core 3.3 supported by wazowski's drivers
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
  #else
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); // I would like this to be 4.1 eventually, but wazowski only has 3.3
  #endif
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
}

template<typename API_t>
void
Application<GLFW_Interface<API_t>>::initialize()
{
  avro_assert( window_.size()!=0 );
  window_[0]->make_current();

  // load gl stuff and print info
  gladLoadGL();
  dumpGLInfo();

  manager_.create_shaders();
}

template<typename API_t>
void
Application<GLFW_Interface<API_t>>::run( const std::string& view )
{
  if (!restart_)
  {
    for (index_t k=0;k<window_.size();k++)
      window_[k]->setup();
  }
  else
    printf("restarting..\n");

  restart_ = false;

  focus_scenes();
  write();

  for (index_t k=0;k<window_.size();k++)
  {
    window_[k]->write_axes();
    window_[k]->clip_plane().set_coordinates( bounding_box_ );
    window_[k]->update_view();
  }

  // option to load existing view
  if (!view.empty())
    window_[0]->controls().load(view);

  // draw everything before entering the event loop
  for (index_t k = 0; k < window_.size(); k++) {
    window_[k]->draw();
    window_[k]->draw(); // draw twice fo ImGui
  }

   // start the event loop
   bool done = false;
   while (!done)
   {
     for (index_t k=0;k<window_.size();k++)
     {
       window_[k]->poll(); // poll for events
       if (window_[k]->should_close())
       {
         printf("window %lu requested close\n",k);
         done = true;
       }
     }

     if (restart_) break;

     #ifdef AVRO_HEADLESS_GRAPHICS
     break;
     #endif
   }

   printf("done\n");

   if (restart_) run();
}

Visualizer::Visualizer()
{
  main_ = std::make_shared<GLFW_Window>(manager_,1024,1024,"avro 2.0 2021");
  add_window( main_.get() );
  //add_window( &side_ );

  initialize();

  main_->create_interface();

  toolbar_ = std::make_shared<Toolbar>(*main_.get(),*this);
  main_->interface().add_widget( toolbar_ );

}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

#endif // AVRO_WITH_GL

} // graphics

} // avro
