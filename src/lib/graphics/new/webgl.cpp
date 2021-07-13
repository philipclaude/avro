
#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"
#include "graphics/new/webgl.h"

namespace avro
{

namespace graphics
{

void
WebGLpp_Manager::add( const VertexAttributeObject& vao ) {

  WebGLpp& gl = webglpp_; // I like to write it like in JavaScript :)

  index_t idx = current_vao_index_;

  // webglpp cannot draw high-order meshes (only high-order solutions)
  avro_assert( vao.order() == 1 );

  // buffer the points
  const PointPrimitive& points = vao.points();
  int vertex_buffer = gl.createBuffer();
  gl.bindBuffer( gl::ARRAY_BUFFER , vertex_buffer );
  gl.bufferData( gl::ARRAY_BUFFER , points.coordinates().data() , sizeof(gl_float) * points.coordinates().size() );
  gl.tagBuffer( gl::ARRAY_BUFFER , "coordinates-" + std::to_string(idx) );

  // buffer the triangles
  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    const TrianglePrimitive& triangles = vao.triangles(k);

    int triangle_buffer = gl.createBuffer();
    printf("triangle_buffer = %d\n",triangle_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , triangle_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , triangles.indices().data() , sizeof(gl_index) * triangles.indices().size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "triangles" + std::to_string(k) + "-" + std::to_string(idx) );
  }

  // buffer the edges
  for (index_t k = 0; k < vao.nb_edges(); k++) {

    const EdgePrimitive& edges = vao.edges(k);

    int edge_buffer = gl.createBuffer();
    printf("edge_buffer = %d\n",edge_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , edge_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , edges.indices().data() , sizeof(gl_index) * edges.indices().size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "edges" + std::to_string(k) + "-" + std::to_string(idx) );
  }

  current_vao_index_++;
}

} // graphics

} // avro
