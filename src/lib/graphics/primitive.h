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

protected:
  coord_t number_;
  const TopologyHolder& topology_;

  int active_;
  ShaderProgram* shader_;
  Window* window_;
  bool visible_;
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
  typedef struct
  {
    std::vector<GLfloat> coordinates;
    std::vector<GLuint> indices;
    std::vector<GLfloat> colours;
    std::vector<GLfloat> normals;
  } glData;

public:
  using Primitive::Primitive;

  void write();
  void draw();
  void convert( glData& data );

private:
  std::vector<GLuint> vbo_;
  GLuint vao_;

  glData data_;

};

} // graphics

} // ursa

#endif
