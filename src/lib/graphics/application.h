#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/types.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace avro
{

class Mesh;
class Primitive;
class Window;

namespace graphics
{

class GraphicsManager
{
public:
  virtual ~GraphicsManager() {}
  virtual void write_points() = 0;
  virtual void write_edges() = 0;
  virtual void write_triangles() = 0;

protected:
};

class OpenGL_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  OpenGL_Manager() {}

  void write_points() { avro_implement; }
  void write_edges() { avro_implement; }
  void write_triangles() { avro_implement; }
};

class Vulkan_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  Vulkan_Manager() {}

  void write_points() { avro_implement; }
  void write_edges() { avro_implement; }
  void write_triangles() { avro_implement; }
};

class WV_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  WV_Manager();

  void write_points() { avro_implement; }
  void write_edges() { avro_implement; }
  void write_triangles() { avro_implement; }
};

class SceneGraph
{
public:
  typedef std::shared_ptr<Primitive> Primitive_ptr;

  void write( GraphicsManager& manager )
  {
    std::map<std::string,Primitive_ptr>::iterator it;
    //for (it=primitive_.begin();it!=primitive_.end();it++)
    //  it->second->write(manager);
  }

private:
  std::map<std::string,Primitive_ptr> primitive_;

  // store all the matrices here
};

class ApplicationBase
{
protected:
  ApplicationBase( GraphicsManager& manager ) :
    manager_(manager)
  {}

  virtual ~ApplicationBase() {}

  virtual void run() = 0;
  void write();

  void add( const Mesh& Mesh );

private:
  GraphicsManager& manager_;
  std::vector<SceneGraph> scene_;
};

template<typename type> class Application;

template<typename T> struct GLFW_Interface;
struct Web_Interface;

template<>
template<typename API_t>
class Application<GLFW_Interface<API_t>> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {}

private:
  API_t manager_;
  std::map<SceneGraph,std::shared_ptr<Window>> scene2window_;
};

template<>
class Application<Web_Interface> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {}

private:
  WV_Manager manager_;
};

} // graphics

} // avro

#endif
