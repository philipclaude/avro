#ifndef avro_LIB_GRAPHICS_PRIMITIVE_H_
#define avro_LIB_GRAPHICS_PRIMITIVE_H_

#include "common/tree.h"

#include "graphics/gl.h"

#include <map>

namespace avro
{

template<typename Master_t> class Topology;
class TopologyBase;

namespace graphics
{

class ShaderProgram;
class Plotter;
class Window;
class GraphicsManager;
class SceneGraph;

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyBase& topology , Window* window );
  Primitive( const TopologyBase& topology , SceneGraph* scene );
  virtual ~Primitive() {}

  void write( GraphicsManager& manager );
  virtual void write() {};
  virtual void draw() {};

  void extract();

  void draw(ShaderProgram&);

  void convert();

  void selectShader( Plotter* plotter );
  ShaderProgram& shader();

  bool& visible() { return visible_; }

  bool& points_on() { return points_on_; }
  bool& edges_on() { return edges_on_; }
  bool& triangles_on() { return triangles_on_; }

  index_t nb_triangles() const { return triangles_.size()/3; }
  index_t nb_edges() const { return edges_.size()/2; }
  index_t nb_points() const { return points_.size()/3; }

  const std::vector<index_t> triangles() const { return triangles0_; }
  const std::vector<index_t> edges() const { return edges0_; }
  const std::vector<real_t> points() const { return points0_; }
  const std::vector<real_t> colors() const { return colors0_; }
  const std::vector<real_t> normals() const { return normals0_; }

  void set_transform_feedback( bool x ) { transform_feedback_ = x; }
  void set_active( const std::string& x ) { active_ = x; }

  coord_t number() const { return number_; }

  Window* window() { return window_; }

protected:
  coord_t number_;
  const TopologyBase& topology_;

  index_t rank_;
  std::string active_;
  ShaderProgram* shader_;
  Window* window_;
  SceneGraph* scene_;
  bool visible_;

  std::vector<GLuint>  edges_;
  std::vector<GLuint>  triangles_;
  std::vector<GLfloat> points_;
  std::vector<GLfloat> normals_;
  std::vector<GLfloat> colors_;

  std::vector<index_t> edges0_;
  std::vector<index_t> triangles0_;
  std::vector<real_t>  points0_;
  std::vector<real_t>  normals0_;
  std::vector<real_t>  colors0_;

  GLuint vao_triangles_;
  GLuint vao_edges_;
  GLuint vao_points_;
  GLuint vao_feedback_;

  bool triangles_on_;
  bool edges_on_;
  bool points_on_;

  bool transform_feedback_;
};

class OpenGLPrimitive : public Primitive
{
private:

public:
  using Primitive::Primitive;

  void draw(GraphicsManager&);
  void write();
  void draw();

};

} // graphics

} // avro

#endif
