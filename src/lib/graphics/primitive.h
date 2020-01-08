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

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyBase& topology , Window* window );
  virtual ~Primitive() {}

  virtual void write() = 0;
  virtual void draw() = 0;

  void selectShader( Plotter* plotter );
  ShaderProgram& shader();

  bool& visible() { return visible_; }

  bool& points_on() { return points_on_; }
  bool& edges_on() { return edges_on_; }
  bool& triangles_on() { return triangles_on_; }

  void set_transform_feedback( bool x ) { transform_feedback_ = x; }
  void set_active( const std::string& x ) { active_ = x; }

protected:
  coord_t number_;
  const TopologyBase& topology_;

  index_t rank_;
  std::string active_;
  ShaderProgram* shader_;
  Window* window_;
  bool visible_;

  std::vector<GLuint>  edges_;
  std::vector<GLuint>  triangles_;
  std::vector<GLfloat> points_;
  std::vector<GLfloat> normals_;
  std::vector<GLfloat> colors_;

  GLuint vao_triangles_;
  GLuint vao_edges_;
  GLuint vao_points_;
  GLuint vao_feedback_;

  bool triangles_on_;
  bool edges_on_;
  bool points_on_;

  bool transform_feedback_;
};

class WebGLPrimitive : public Primitive
{
public:
  using Primitive::Primitive;

  void write();
  void draw();

private:
  int handle_;

};

class OpenGLPrimitive : public Primitive
{
private:

public:
  using Primitive::Primitive;

  void write();
  void draw();
  void convert();

};

} // graphics

} // avro

#endif
