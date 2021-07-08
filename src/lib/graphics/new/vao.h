#ifndef AVRO_LIB_GRAPHICS_VAO_H_
#define AVRO_LIB_GRAPHICS_VAO_H_

#include "avro_types.h"

#include "common/error.h"

#include "graphics/gl.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include <vector>

namespace avro
{

namespace graphics
{

class PointPrimitive {

public:

	PointPrimitive( const Points& points ) :
    points_(points) {

    coord_t nb_dim = points.dim();
    if (nb_dim > 3) nb_dim = 3;
    avro_assert( points.nb_ghost() == 0 );
    for (index_t k = 0; k < points.nb(); k++) {
      for (coord_t d = 0; d < nb_dim; d++)
        coordinates_.push_back( points[k][d] );
      for (coord_t d = nb_dim; d < 3; d++)
        coordinates_.push_back( 0.0 );
    }
  }

  ~PointPrimitive() {
    glDeleteBuffers( 1 , &buffer_ );
  }

	void write() {
    GL_CALL( glGenBuffers( 1 , &buffer_ ) );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * coordinates_.size() , coordinates_.data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );


  }

  void draw( ShaderProgram& program ) {
    program.use();
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , buffer_) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , coordinates_.size()/3 ) );
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
  }

	gl_index buffer() const { return buffer_; }

  index_t nb() const { return coordinates_.size() / 3; }

  void print() const {
    for (index_t k = 0; k < nb(); k++) {
      printf("point %lu = ( ",k);
      for (coord_t d = 0; d < 3; d++)
        printf("%g ",coordinates_[3*k+d]);
      printf(")\n");
    }
  }

private:
  const Points& points_;
	std::vector<gl_float> coordinates_;
	gl_index buffer_;
	bool  visible_;
};

class EdgePrimitive {

public:
  EdgePrimitive( coord_t order ) :
    order_(order),
    nb_basis_(order+1)
  {}

  void add( const index_t* v , index_t nv ) {
    avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
    for (index_t j = 0; j < nv; j++)
      indices_.push_back(v[j]);
  }

  index_t nb() const { return indices_.size() / nb_basis_; }

  void write() {
    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &buffer_ ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * indices_.size() , indices_.data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
  }

  void draw() {
    // bind the buffer for the indices we want to draw
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }

  void print() const {
    for (index_t k = 0; k < nb(); k++) {
      printf("edge %lu = ( ",k);
      for (index_t j = 0; j < nb_basis_; j++)
        printf("%d ",int(indices_[nb_basis_*k+j]));
      printf(")\n");
    }
  }

private:
  coord_t order_;
  index_t nb_basis_;
  gl_index buffer_;
  std::vector<gl_index> indices_;
};

class FieldPrimitive;

class TrianglePrimitive {
public:
  TrianglePrimitive( coord_t order ) :
    order_(order),
    nb_basis_(nb_simplex_basis(2,order_)),
    field_(nullptr)
  {}

  ~TrianglePrimitive() {
    glDeleteBuffers( 1 , &buffer_ );
  }

  void add( const index_t* v , index_t nv ) {
    avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
    for (index_t j = 0; j < nv; j++)
      indices_.push_back(v[j]);
  }

  void print() const {
    for (index_t k = 0; k < nb(); k++) {
      printf("triangle %lu = ( ",k);
      for (index_t j = 0; j < nb_basis_; j++)
        printf("%d ",int(indices_[nb_basis_*k+j]));
      printf(")\n");
    }
  }

  index_t nb() const { return indices_.size() / nb_basis_; }

  void write() {
    // bind the triangles
    GL_CALL( glGenBuffers( 1 , &buffer_ ) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * indices_.size() , indices_.data() , GL_STATIC_DRAW) );
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
  }

  void draw();

  void set_field( FieldPrimitive* field ) { field_ = field; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_index> indices_;
  gl_index buffer_;
  FieldPrimitive* field_;
};

class FieldData {
public:
  FieldData( coord_t order ) :
    order_(order),
    nb_basis_(nb_simplex_basis(2,order_))
  {}

  void add( real_t* f , index_t ndof ) {
    avro_assert( ndof == nb_basis_ );
    for (index_t j = 0; j < ndof; j++)
      data_.push_back(f[j]);
  }

  const std::vector<gl_float>& data() const { return data_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_float> data_;
};

class FieldPrimitive {
public:
  FieldPrimitive() {}

  ~FieldPrimitive() {

  }

  void set_active( const std::string& active ) {
    active_ = active;
  }

  void add( const std::string& name , std::shared_ptr<FieldData> data ) {
    data_.insert({name,data});
  }

  void write() {

    FieldData* field = data_[active_].get();

    GL_CALL( glGenBuffers( 1 , &buffer_ ) );
    GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , buffer_ ) );
    GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(gl_float) * field->data().size() , field->data().data() , GL_STATIC_DRAW) );

    GL_CALL( glGenTextures( 1 , &texture_) );
    GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) ); // fields are always in texture 0
    GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture_) );
    GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , buffer_ ) );
  }

  void activate( ShaderProgram& shader ) {
    // bind the desired solution texture
    glActiveTexture(GL_TEXTURE0 + 0);
    GLint solution_location = glGetUniformLocation(shader.handle() , "solution");
    //avro_assert( solution_location >= 0 );
    glUniform1i(solution_location, 0); // first sampler in fragment shader
  }

private:
  std::string active_;
  std::map<std::string, std::shared_ptr<FieldData> > data_;
  gl_index texture_;
  gl_index buffer_;
};

struct MeshFacet;

class VertexAttributeObject {

public:
  VertexAttributeObject( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  void build( const TopologyBase& topology );

  void draw_triangles( ShaderProgram& );
  void draw_edges( ShaderProgram& );
  void draw_points( ShaderProgram& );

private:
  template<typename type>
  void get_primitives( const Topology<type>& topology , const std::vector<std::vector<MeshFacet>>& facets );

  template<typename type>
  void _build( const Topology<type>& topology );

private:
  coord_t number_;
  coord_t order_;

  gl_index vertex_array_;
  std::shared_ptr<PointPrimitive> points_;
  std::vector< std::shared_ptr<EdgePrimitive> > edges_;
  std::vector< std::shared_ptr<TrianglePrimitive> > triangles_;
  std::vector< std::shared_ptr<FieldPrimitive> > solution_;
};

} // graphics

} // avro

#endif
