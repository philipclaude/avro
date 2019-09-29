#ifndef URSA_LIB_GRAPHICS_PRIMITIVE_H_
#define URSA_LIB_GRAPHICS_PRIMITIVE_H_

#include "common/tree.h"

namespace ursa
{

template<typename Master_t> class Topology;
class TopologyHolder;
class Fields;

namespace graphics
{

class Primitive : public Tree<Primitive>
{
public:
  Primitive( const TopologyHolder& topology , const Fields* fields=nil );
  virtual ~Primitive() {}

  virtual void write() = 0;

protected:
  const TopologyHolder& topology_;
  const Fields* fields_;

  int identifier_;
  int active_;
};

class WebGLPrimitive : public Primitive
{
public:
  using Primitive::Primitive;

  void write();

};

class OpenGLPrimitive : public Primitive
{
public:
  using Primitive::Primitive;

  void write();

private:

};

} // graphics

} // ursa

#endif
