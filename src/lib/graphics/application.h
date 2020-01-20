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

  void set_matrices( SceneGraph& scene );

private:

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
};

class TransformFeedbackResult
{
public:
  GLuint& vao() { return vao_; }
  GLuint& query() { return query_; }
  GLuint& buffer() { return buffer_; }

private:
  GLuint vao_;
  GLuint query_;
  GLuint buffer_;

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

  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr );

private:

  void set_matrices( SceneGraph& scene );
  void draw( Primitive& primitive , TransformFeedbackResult* feedback=nullptr );

  // map from Primitive to vao
  std::map<Primitive*,index_t> vao_points_;
  std::map<Primitive*,index_t> vao_edges_;
  std::map<Primitive*,index_t> vao_triangles_;

  std::map<Primitive*,ShaderProgram*> shader_;
  ShaderLibrary shaders_;
};

class Vulkan_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  Vulkan_Manager() {}

  void write_points() { avro_implement; }
  void write_edges() { avro_implement; }
  void write_triangles() { avro_implement; }

  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_implement; }

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

  SceneGraph() :
    update_(true)
  {}

  void add_primitive( TopologyBase& topology , const std::string& shader_name="default" )
  {
    Primitive_ptr primitive ;avro_implement;//= std::make_shared<Primitive>(topology,shader,*this);
    primitive_.push_back(primitive);
  }

  void write( GraphicsManager& manager )
  {
    for (index_t k=0;k<primitive_.size();k++)
      primitive_[k]->write(manager);
  }

  bool update() const { return update_; }
  void set_update( bool x ) { update_ = x; }

  const mat4& mvp_matrix() const { return mvpMatrix_; }
  const mat4& normal_matrix() const { return normalMatrix_; }

  mat4& mvp_matrix() { return mvpMatrix_; }
  mat4& normal_matrix() { return normalMatrix_; }

  index_t nb_primitives() const { return primitive_.size(); }

  Primitive& primitive( index_t k ) { return *primitive_[k].get(); }

private:
  std::vector<Primitive_ptr> primitive_; // roots of the scene graph

  // store all the matrices here
  mat4 mvpMatrix_;
  mat4 viewMatrix_;
  mat4 projMatrix_;
  mat4 modelMatrix_;
  mat4 normalMatrix_;
  mat4 modelViewMatrix_;

  bool update_;
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
    scene_.push_back(SceneGraph());
  }

  SceneGraph& scene( index_t k ) { return scene_[k]; }

  index_t nb_scene() const { return scene_.size(); }

private:
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

  void run()
  {
     // start the rendering loop
     while (true)
     {
       for (index_t k=0;k<nb_scene();k++)
        manager_.draw(scene(k));
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
    manager_gl.draw(scene(0),&feedback);
  }

private:
  WV_Manager manager_;
};

} // graphics

} // avro

#endif
