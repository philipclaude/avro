//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/primitives.h"
#include "graphics/vao.h"
#include "graphics/webglpp.h"

#include <iostream>
#include <vector>
#include <math.h>

namespace avro
{

namespace graphics
{

void
add_json_field( std::string& J , const std::string& name , const int& i , bool end=false) {
  J += "\u0022" + name + "\u0022:" + std::to_string(i);
  if (end) J += "}";
  else J += ",";
}

void
add_json_field( std::string& J , const std::string& name , const std::string& s , bool end=false) {
  J += "\u0022" + name + "\u0022:\u0022" + s + "\u0022";
  if (end) J += "}";
  else J += ",";
}

template<typename T>
void
add_json_field( std::string& J , const std::string& name , const std::vector<T>& x , bool end=false) {
  J += "\u0022" + name + "\u0022:[";
  if (x.size() == 0)
    J += "]";
  else {
    for (int k = 0; k < (int)x.size(); k++) {
      J += std::to_string(x[k]);
      if (k < int(x.size())-1) J += ",";
      else J += "]";
    }
  }
  if (end) J += "}";
  else J += ",";
}

void
WebGLpp::send( int port ) {

  std::string message = "{ \"buffers\": [";
  for (int i = 0; i < (int)buffers_.size(); i++) {

    const BufferObject& buffer = *buffers_[i].get();

    // retrieve the data
    std::string J = "{";
    if (buffer.type() == typeid(float).name()) {
      std::vector<float> x = buffer.Float32Array();
      add_json_field(J,"data",x);
      add_json_field(J,"type","Float32Array");
    }
    else if (buffer.type() == typeid(unsigned short).name()) {
      std::vector<unsigned short> x = buffer.Uint16Array();
      add_json_field(J,"data",x);
      add_json_field(J,"type","Uint16Array");
    }
    else if (buffer.type() == typeid(unsigned int).name()) {
      std::vector<unsigned int> x = buffer.Uint32Array();
      add_json_field(J,"data",x);
      add_json_field(J,"type","Uint32Array");
    }
    else {
      ERROR("unimplemented buffer type");
    }
    add_json_field(J,"tag",buffer.tag());
    add_json_field(J,"index",i);
    add_json_field(J,"target",target_name(buffer.target()),true);

    message += J;
    if (i < (int)buffers_.size()-1) message += ",";
    else message += "]}";
  }

  // at the end we write the data to the client
  websockets::write( port , {message} );
}

void
WebGL_Manager::write( const VertexArrayObject& vao ) {

  WebGLpp& gl = webglpp_; // I like to write it like in JavaScript :)

  index_t idx = current_vao_index_;

  // webglpp cannot draw high-order meshes (only high-order solutions)
  avro_assert( vao.order() == 1 );

  if (vao.nb_fields() == 0) {

    // buffer the triangles
    for (index_t k = 0; k < vao.nb_triangles(); k++) {

      const TrianglePrimitive& triangles = vao.triangles(k);

      int triangle_buffer = gl.createBuffer();
      gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , triangle_buffer );
      gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , triangles.indices().data() , sizeof(gl_index) * triangles.indices().size() );
      gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "triangles" + std::to_string(k) + "-vao" + std::to_string(idx) );
    }

    // buffer the points
    const PointPrimitive& points = vao.points();
    int vertex_buffer = gl.createBuffer();
    gl.bindBuffer( gl::ARRAY_BUFFER , vertex_buffer );
    gl.bufferData( gl::ARRAY_BUFFER , points.coordinates().data() , sizeof(gl_float) * points.coordinates().size() );
    gl.tagBuffer( gl::ARRAY_BUFFER , "coordinates-vao" + std::to_string(idx) );


    // buffer the edges
    for (index_t k = 0; k < vao.nb_edges(); k++) {

      const EdgePrimitive& edges = vao.edges(k);

      int edge_buffer = gl.createBuffer();
      gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , edge_buffer );
      gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , edges.indices().data() , sizeof(gl_index) * edges.indices().size() );
      gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "edges" + std::to_string(k) + "-vao" + std::to_string(idx) );
    }

  }
  else {
    // we need to duplicate points for the field rendering to work
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
          if (point_map.find(vertex) == point_map.end())
            point_map.insert( {vertex,n} );
          indices.push_back(n++);
        }
      }

      int triangle_buffer = gl.createBuffer();
      gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , triangle_buffer );
      gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , indices.data() , sizeof(gl_index) * indices.size() );
      gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "triangles" + std::to_string(k) + "-vao" + std::to_string(idx) );
    }

    // buffer the points
    int vertex_buffer = gl.createBuffer();
    gl.bindBuffer( gl::ARRAY_BUFFER , vertex_buffer );
    gl.bufferData( gl::ARRAY_BUFFER , coordinates.data() , sizeof(gl_float) * coordinates.size() );
    gl.tagBuffer( gl::ARRAY_BUFFER , "coordinates-vao" + std::to_string(idx) );

    // buffer the edges
    for (index_t k = 0; k < vao.nb_edges(); k++) {

      const std::vector<gl_index>& edges = vao.edges(k).indices();

      std::vector<gl_index> indices( edges.size() );
      for (index_t j = 0; j < edges.size(); j++)
        indices[j] = point_map.at(edges[j]);

      int edge_buffer = gl.createBuffer();
      gl.bindBuffer( gl::ELEMENT_ARRAY_BUFFER , edge_buffer );
      gl.bufferData( gl::ELEMENT_ARRAY_BUFFER , indices.data() , sizeof(gl_index) * indices.size() );
      gl.tagBuffer( gl::ELEMENT_ARRAY_BUFFER , "edges" + std::to_string(k) + "-vao" + std::to_string(idx) );
    }

    // buffer the fields
    for (index_t k = 0; k < vao.nb_fields(); k++) {

      const FieldPrimitive& solution = vao.field(k);
      const std::map<std::string,std::shared_ptr<FieldData>>& data = solution.data();
      std::map<std::string,std::shared_ptr<FieldData>>::const_iterator it; // for now just retrieve field in first rank

      for (it = data.begin(); it != data.end(); ++it) {
        const FieldData& field = *it->second.get();
        const std::string& name = it->first;

        int field_buffer = gl.createBuffer();
        gl.bindBuffer( gl::ARRAY_BUFFER , field_buffer );
        gl.bufferData( gl::ARRAY_BUFFER , field.data().data() , sizeof(gl_float) * field.data().size() );
        gl.tagBuffer( gl::ARRAY_BUFFER , "field_" + name + "_order=" + std::to_string(field.order()) + "_tri" + std::to_string(k) + "-vao" + std::to_string(idx) );
      }
    }
  }

  current_vao_index_++;
}


} // graphics

} // avro
