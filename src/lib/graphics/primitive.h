#ifndef URSA_LIB_GRAPHICS_PRIMITIVE_H_
#define URSA_LIB_GRAPHICS_PRIMITIVE_H_

#include "common/tree.h"

#include "graphics/gl.h"

#include <map>

namespace ursa
{

template<typename Master_t> class Topology;
class TopologyHolder;

namespace graphics
{

class ShaderProgram;
class Plotter;
class Window;

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyHolder& topology , Window* window );
  virtual ~Primitive() {}

  virtual void write() = 0;
  virtual void draw() = 0;

  void selectShader( Plotter* plotter );
  ShaderProgram& shader();

  bool& visible() { return visible_; }

  bool& points_on() { return points_on_; }
  bool& edges_on() { return edges_on_; }
  bool& triangles_on() { return triangles_on_; }

protected:
  coord_t number_;
  const TopologyHolder& topology_;

  int active_;
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

  bool triangles_on_;
  bool edges_on_;
  bool points_on_;
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
  /*typedef struct
  {
    std::vector<GLfloat> coordinates;
    std::vector<GLuint> indices;
    std::vector<GLfloat> colours;
    std::vector<GLfloat> normals;
  } glData;*/

public:
  using Primitive::Primitive;

  void write();
  void draw();
  void convert();

private:
  //std::vector<GLuint> vbo_;
  //GLuint vao_;

  //glData data_;

};

} // graphics

} // ursa

#endif
