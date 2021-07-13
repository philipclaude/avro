#ifndef AVRO_LIB_GRAPHICS_WEBGL_H_
#define AVRO_LIB_GRAPHICS_WEBGL_H_

#include "avro_types.h"

#include "webglpp.h"

namespace avro
{

namespace graphics
{

class VertexAttributeObject;

class WebGLpp_Manager {

public:
  WebGLpp_Manager( int port ) :
    webglpp_(port),
    current_vao_index_(0)
  {}

  void add( const VertexAttributeObject& vao );

  void send() {
    webglpp_.send();
  }

private:
  WebGLpp webglpp_;
  index_t current_vao_index_;
};

}

} // avro

#endif
