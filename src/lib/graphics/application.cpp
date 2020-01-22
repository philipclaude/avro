#include "graphics/application.h"
#include "graphics/gl.h"
#include "graphics/window.h"

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
  for (index_t k=0;k<window_.size();k++)
    window_[k]->setup();

  write();

  for (index_t k=0;k<window_.size();k++)
    window_[k]->update_view();

  const double fps = 1.0 / 60.0;
  double last_update_time = 0;  // number of seconds since the last loop
  double last_frame_time = 0;   // number of seconds since the last frame

   // start the rendering loop
   bool done = false;
   while (!done)
   {

     double now = glfwGetTime();

     // This if-statement only executes once every 60th of a second
     if ((now - last_frame_time) >= fps)
     {
       // draw your frame here
       for (index_t k=0;k<window_.size();k++)
       {
         window_[k]->begin_draw();

         for (index_t j=0;j<window_[k]->nb_scene();j++)
           manager_.draw(window_[k]->scene(j));

         window_[k]->end_draw();

         if (window_[k]->should_close())
         {
           printf("window %lu requested close\n",k);
           done = true;
         }

       }

       // only set lastFrameTime when you actually draw something
       last_frame_time = now;
      }

      // set lastUpdateTime every iteration
      last_update_time = now;
   }
}

Visualizer::Visualizer()
{
  main_ = std::make_shared<GLFW_Window>(manager_,1024,612,"3D");
  add_window( main_.get() );
  //add_window( &side_ );

  initialize();
}

void
Visualizer::add_topology( const TopologyBase& topology )
{
  // add the topology to the relevant windows
  index_t id = main_->create_scene();
  scenes_.push_back( &main_->scene(id) );
  index_t prim_id = main_->scene(id).add_primitive(topology); // create a new root in the scene graph
  manager_.select_shader( main_->scene(id).primitive(prim_id) , "wv" );
}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

} // graphics

} // avro
