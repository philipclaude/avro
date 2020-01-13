#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include "graphics/primitive.h"
#include "graphics/shader.h"

#include <map>
#include <memory>
#include <unordered_set>
#include <string>
#include <vector>

namespace avro
{

class Mesh;
class Window;

namespace graphics
{

class SceneGraph;

class EmptyShader : public ShaderProgram
{
  // used for webviewer
};

class ShaderLibrary
{
public:
  ShaderLibrary()
  {
    // generate all the shaders!!
  }

  const ShaderProgram& operator[] ( const std::string& name ) const
  {
    avro_assert( shaders_.find(name)!=shaders_.end() );
    return shaders_.at(name);
  }

  ShaderProgram& operator[] ( const std::string& name )
  {
    avro_assert( shaders_.find(name)!=shaders_.end() );
    return shaders_.at(name);
  }

  // store all the shaders
  std::map<std::string,ShaderProgram> shaders_;
};

class GraphicsManager
{
public:
  virtual ~GraphicsManager() {}
  virtual void write_points() = 0;
  virtual void write_edges() = 0;
  virtual void write_triangles() = 0;

protected:
  ShaderLibrary shaders_;
};

class TransformFeedbackResult
{

};

class OpenGL_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  OpenGL_Manager()
  {}

public:
  void write_points() { avro_implement; }
  void write_edges() { avro_implement; }
  void write_triangles() { avro_implement; }

  void draw();
  void capture( SceneGraph& scene , TransformFeedbackResult& feedback );

private:

  // map from vao to number of primitives to draw
  std::map<index_t,index_t> vao_points_;
  std::map<index_t,index_t> vao_edges_;
  std::map<index_t,index_t> vao_triangles_;

  // map from vao to whether the primitive is active
  std::map<index_t,bool> vao_points_active_;
  std::map<index_t,bool> vao_edges_active_;
  std::map<index_t,bool> vao_triangles_active_;

  // map from vao to shader program
  std::map<index_t,ShaderProgram*> vao_points_shader_;
  std::map<index_t,ShaderProgram*> vao_edges_shader_;
  std::map<index_t,ShaderProgram*> vao_triangles_shader_;
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

  SceneGraph( ShaderLibrary& shaders ) :
    shader_library_(shaders),
    needs_update_(false)
  {}

  void add_primitive( TopologyBase& topology , const std::string& shader_name="default" )
  {
    ShaderProgram* shader = &shader_library_[shader_name];
    Primitive_ptr primitive ;avro_implement;//= std::make_shared<Primitive>(topology,shader,*this);
    primitive_.push_back(primitive);
    shaders_.push_back(shader);
    uniquify(shaders_);
  }

  void write( GraphicsManager& manager )
  {
    for (index_t k=0;k<primitive_.size();k++)
      primitive_[k]->write(manager);
  }

  void set_matrices()
  {
    // go through all the active shaders and assign the MVP and normalMatrix
    for (index_t k=0;k<shaders_.size();k++)
    {
      ShaderProgram& shader = *shaders_[k];
      shader.setUniform("MVP" , mvpMatrix_ );
      shader.setUniform("u_normalMatrix" , normalMatrix_ );
    }
  }

  void draw()
  {
    if (!needs_update_) return;

    set_matrices();

    for (index_t k=0;k<primitive_.size();k++)
    {
      // use the shader to draw
      primitive_[k]->draw();
    }

    needs_update_ = false;
  }

private:
  std::vector<Primitive_ptr> primitive_;
  std::vector<ShaderProgram*> shaders_; // active shaders <= primitive_.size()
  ShaderLibrary& shader_library_; // list of all shaders

  // store all the matrices here
  mat4 mvpMatrix_;
  mat4 viewMatrix_;
  mat4 projMatrix_;
  mat4 modelMatrix_;
  mat4 normalMatrix_;
  mat4 modelViewMatrix_;

  bool needs_update_;
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

  void add_scene()
  {
    scene_.push_back(SceneGraph(shader_library_));
  }
  SceneGraph& scene( index_t k ) { return scene_[k]; }

  index_t nb_scene() const { return scene_.size(); }

private:
  ShaderLibrary shader_library_;
  GraphicsManager& manager_;
  std::vector<SceneGraph> scene_;
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
  {}

  void render();

  void run()
  {
     // start the rendering loop
     while (true)
     {
       for (index_t k=0;k<nb_scene();k++)
        scene(k).draw();
     }
  }

private:
  API_t manager_;
  std::vector<std::shared_ptr<Window>> window_;
};

template<>
class Application<Web_Interface> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_)
  {
    add_scene(); // single scene graph passed to the web for now
  }

  void run() { avro_implement; } // run the server

  void save_eps() // will also need transformation matrices sent from the web
  {
    // first set the transformation to the scene
    OpenGL_Manager manager_gl;
    scene(0).write(manager_gl);
    TransformFeedbackResult feedback;
    manager_gl.capture(scene(0),feedback);
  }

private:
  WV_Manager manager_;
};

} // graphics

} // avro

#endif
