#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include "graphics/gl.h"
#include "graphics/manager.h"
#include "graphics/scene.h"
#include "graphics/window.h"
#include "graphics/wv.h"

#include <map>
#include <memory>
#include <unordered_set>
#include <string>
#include <vector>

namespace avro
{

class Mesh;
class TopologyBase;

namespace graphics
{

class ApplicationBase
{
public:
  void receive( const std::string& text ) const;

  void write();

protected:
  ApplicationBase( GraphicsManager& manager ) :
    manager_(manager)
  {}

  virtual ~ApplicationBase() {}

  virtual void run() = 0;

protected:
  std::vector<SceneGraph*> scenes_;

private:
  GraphicsManager& manager_;
};

template<typename type> class Application;

template<typename T> struct GLFW_Interface;
struct Web_Interface;

template<typename API_t>
class Application<GLFW_Interface<API_t>> : public ApplicationBase
{
public:
  Application() :
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
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    #endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  }

  void initialize();
  void run();

  bool& restart() { return restart_; }

protected:
  void add_window( GLFW_Window* window )
    { window_.push_back(window); }

protected:
  API_t manager_;
  std::vector<GLFW_Window*> window_;
  bool restart_;
};

template<>
class Application<Web_Interface> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {}

  void run();
  void save_eps();

protected:
  WV_Manager manager_;
  SceneGraph scene_;

private:
  void connect_client();
};

class Visualizer : public Application<GLFW_Interface<OpenGL_Manager>>
{
public:
  Visualizer();

  void add_topology( const TopologyBase& topology )
  {
    //static_cast< const Topology<Simplex>& >(topology).Tree<Topology<Simplex>>::print();

    // add the topology to the relevant windows
    index_t id = main_->create_scene();
    scenes_.push_back( &main_->scene(id) );
    index_t prim_id = main_->scene(id).add_primitive(topology); // create a new root in the scene graph
    manager_.select_shader( main_->scene(id).primitive(prim_id) , "wv" );

    for (index_t j=0;j<main_->scene(id).primitive(prim_id).nb_children();j++)
      manager_.select_shader( main_->scene(id).primitive(prim_id).child(j) , "wv" );
  }

  void remove( index_t scene , index_t root )
  {
    scenes_[scene]->remove(root);
    if (scenes_[scene]->nb_primitives()==0)
    {
      scenes_.erase( scenes_.begin() + root );
    }
  }

  OpenGL_Manager& manager() { return manager_; }

  GLFW_Window& main_window() { return *main_.get(); }

  std::shared_ptr<GLFW_Window> main_;
  //GLFW_Window side_;
};

class WebVisualizer : public Application<Web_Interface>
{
public:

  WebVisualizer()
  {
    scenes_.push_back( &scene_ );
  }

  template<typename type>
  void add_topology( Topology<type>& topology )
  {
    // create a new root in the scene graph
    scene_.add_primitive(topology);
  }

private:
  //SceneGraph scene_;
};

} // graphics

} // avro

#endif
