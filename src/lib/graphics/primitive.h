#ifndef URSA_LIB_GRAPHICS_PRIMITIVE_H_
#define URSA_LIB_GRAPHICS_PRIMITIVE_H_

#include "common/tree.h"

#include "graphics/gl.h"

#include <map>

namespace ursa
{

template<typename Master_t> class Topology;
class TopologyHolder;
class Fields;

namespace graphics
{

class ShaderProgram;
class Plotter;

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyHolder& topology , const Fields* fields=nil );
  virtual ~Primitive() {}

  virtual void write() = 0;
  virtual void draw() = 0;

  void selectShader( Plotter* plotter );
  ShaderProgram& shader();

protected:
  coord_t number_;
  const TopologyHolder& topology_;
  const Fields* fields_;

  int active_;
  ShaderProgram* shader_;
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
public:
  using Primitive::Primitive;

  void write();
  void draw();

private:
  std::vector<GLuint> vbo_;
  GLuint vao_;

};

} // graphics

} // ursa

#endif
