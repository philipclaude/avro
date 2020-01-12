#ifndef avro_LIB_GRAPHICS_SCENE_H_
#define avro_LIB_GRAPHICS_SCENE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace avro
{

namespace graphics
{

class Plot;
class Shader;
class Primitive;

class SceneGraph
{
  typedef std::shared_ptr<Primitive> Primitive_ptr;

public:

  void add_primitive( Primitive_ptr prim )
  {
    primitive_.insert( {"0x1" , prim} );
  }

private:
  // you might be wondering: why store a string representation of
  // the pointer to the primitive? why not just use a vector of
  // primitive and look up a pointer as necessary?
  // well, when sending things over to the WebViewer, we can't
  // compare pointers anymore...we need to look things up via some
  // identifier, hence a string
  std::map<std::string,Primitive_ptr> primitive_;
};

} // graphics

} // avro

#endif
