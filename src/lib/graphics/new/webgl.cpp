
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

  #if 0
  // buffer the triangles
  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    const TrianglePrimitive& triangles = vao.triangles(k);

    int triangle_buffer = gl.createBuffer();
    printf("triangle_buffer = %d\n",triangle_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , triangle_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , triangles.indices().data() , sizeof(gl_index) * triangles.indices().size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "triangles" + std::to_string(k) + "-" + std::to_string(idx) );
  }

  // buffer the points
  const PointPrimitive& points = vao.points();
  int vertex_buffer = gl.createBuffer();
  gl.bindBuffer( gl::ARRAY_BUFFER , vertex_buffer );
  gl.bufferData( gl::ARRAY_BUFFER , points.coordinates().data() , sizeof(gl_float) * points.coordinates().size() );
  gl.tagBuffer( gl::ARRAY_BUFFER , "coordinates-" + std::to_string(idx) );


  // buffer the edges
  for (index_t k = 0; k < vao.nb_edges(); k++) {

    const EdgePrimitive& edges = vao.edges(k);

    int edge_buffer = gl.createBuffer();
    printf("edge_buffer = %d\n",edge_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , edge_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , edges.indices().data() , sizeof(gl_index) * edges.indices().size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "edges" + std::to_string(k) + "-" + std::to_string(idx) );
  }


  #else

  std::map<gl_index,gl_index> point_map;
  const std::vector<gl_float>& points = vao.points().coordinates();
  std::vector<gl_float> coordinates;
  index_t n = 0;
  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    const std::vector<gl_index>& triangles = vao.triangles(k).indices();
    std::vector<gl_index> indices;
    for (index_t j = 0 ; j < triangles.size()/3; j++) {
      for (coord_t i = 0; i < 3; i++) {
        index_t vertex = triangles[3*j+i];
        for (coord_t d = 0; d < 3; d++)
          coordinates.push_back( points[3*vertex + d] );
        point_map.insert( {vertex,n} );
        indices.push_back(n++);
      }
    }

    int triangle_buffer = gl.createBuffer();
    printf("triangle_buffer = %d\n",triangle_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , triangle_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , indices.data() , sizeof(gl_index) * indices.size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "triangles" + std::to_string(k) + "-" + std::to_string(idx) );
  }

  // buffer the points
  int vertex_buffer = gl.createBuffer();
  gl.bindBuffer( gl::ARRAY_BUFFER , vertex_buffer );
  gl.bufferData( gl::ARRAY_BUFFER , coordinates.data() , sizeof(gl_float) * coordinates.size() );
  gl.tagBuffer( gl::ARRAY_BUFFER , "coordinates-" + std::to_string(idx) );

  // buffer the edges
  for (index_t k = 0; k < vao.nb_edges(); k++) {

    const std::vector<gl_index>& edges = vao.edges(k).indices();

    std::vector<gl_index> indices;
    for (index_t j = 0; j < edges.size(); j++)
      indices.push_back( point_map[edges[j]] );

    int edge_buffer = gl.createBuffer();
    printf("edge_buffer = %d\n",edge_buffer);
    gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , edge_buffer );
    gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , indices.data() , sizeof(gl_index) * indices.size() );
    gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "edges" + std::to_string(k) + "-" + std::to_string(idx) );
  }


  #endif

  // buffer the fields
  for (index_t k = 0; k < vao.nb_fields(); k++) {

    const FieldPrimitive& solution = vao.field(k);
    const std::map<std::string,std::shared_ptr<FieldData>>& data = solution.data();
    const std::map<std::string,std::shared_ptr<FieldData>>::const_iterator it = data.begin(); // for now just retrieve field in first rank
    const FieldData& field = *it->second.get();
    //const std::string& name = it->first;

    int field_buffer = gl.createBuffer();
    printf("field_buffer = %d\n",field_buffer);
    gl.bindBuffer( gl::ARRAY_BUFFER , field_buffer );
    gl.bufferData( gl::ARRAY_BUFFER , field.data().data() , sizeof(gl_float) * field.data().size() );
    gl.tagBuffer( gl::ARRAY_BUFFER , "field_order=" + std::to_string(field.order()) + "_tri" + std::to_string(k) + "-" + std::to_string(idx) );
  }

  current_vao_index_++;
}

} // graphics

} // avro
