#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include "graphics/gl.h"
#include "graphics/manager.h"
#include "graphics/scene.h"
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

protected:
  ApplicationBase( GraphicsManager& manager ) :
    manager_(manager)
  {}

  virtual ~ApplicationBase() {}

  virtual void run() = 0;
  void write();

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
    ApplicationBase(manager_)
  {
    // initialize OpenGL
    avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );

    // set the version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  }

  void initialize();
  void run();

protected:
  void add_window( GLFW_Window* window )
    { window_.push_back(window); }

protected:
  API_t manager_;
  std::vector<GLFW_Window*> window_;
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

  void add_topology( const TopologyBase& topology );

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

  void add_topology( TopologyBase& topology )
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
