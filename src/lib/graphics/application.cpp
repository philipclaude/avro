#include "graphics/application.h"
#include "graphics/gl.h"
#include "graphics/window.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

namespace avro
{

namespace graphics
{

void
ApplicationBase::write()
{
  // write all the data present in the scene graph
  printf("nb_scenes = %lu\n",scenes_.size());
  for (index_t k=0;k<scenes_.size();k++)
    scenes_[k]->write(manager_);
}

void
ApplicationBase::receive( const std::string& text ) const
{
  printf("received text %s\n",text.c_str());
}

void
Application<Web_Interface>::save_eps()
{
  // first set the transformation to the scene
  OpenGL_Manager manager_gl;
  scene_.write(manager_gl);
  TransformFeedbackResult feedback;
  manager_gl.draw(scene_,&feedback);
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
Application<GLFW_Interface<API_t>>::run()
{
  if (!restart_)
  {
    for (index_t k=0;k<window_.size();k++)
      window_[k]->setup();
  }
  else
    printf("restarting..\n");

  restart_ = false;

  write();

  for (index_t k=0;k<window_.size();k++)
    window_[k]->update_view();

  index_t fps = 120;
  for (index_t k=0;k<window_.size();k++)
    fps = std::min( fps , window_[k]->fps() );

  const double spf = 1.0 / fps;
  double last_update_time = 0;  // number of seconds since the last loop
  double last_frame_time = 0;   // number of seconds since the last frame
  UNUSED(last_update_time);

   // start the rendering loop
   bool done = false;
   while (!done)
   {

     double now = glfwGetTime();

     // this if-statement only executes once every 60th of a second
     if ((now - last_frame_time) >= spf)
     {
       // draw frame
       for (index_t k=0;k<window_.size();k++)
       {
         window_[k]->begin_draw();
         window_[k]->update_view();

         if (!restart_)
         for (index_t j=0;j<window_[k]->nb_scene();j++)
           manager_.draw(window_[k]->scene(j));

         window_[k]->end_draw();

         if (window_[k]->should_close())
         {
           printf("window %lu requested close\n",k);
           done = true;
         }

       }

       if (restart_) break;

       // only set lastFrameTime when you actually draw something
       last_frame_time = now;
      }

      #ifdef AVRO_HEADLESS_GRAPHICS
      break;
      #endif

      // set lastUpdateTime every iteration
      last_update_time = now;
   }

   if (restart_) run();
}

Visualizer::Visualizer()
{
  main_ = std::make_shared<GLFW_Window>(manager_,1024,1024,"avro 2.0 2020");
  add_window( main_.get() );
  //add_window( &side_ );

  initialize();

  main_->create_interface();

}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

} // graphics

} // avro
