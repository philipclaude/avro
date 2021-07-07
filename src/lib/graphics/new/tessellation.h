#ifndef AVRO_LIB_GRAPHICS_MESH_H_
#define AVRO_LIB_GRAPHICS_MESH_H_

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

  void render() {
		// draw the points
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

class TrianglePrimitive {
public:
  TrianglePrimitive( coord_t order ) :
    order_(order),
    nb_basis_(nb_simplex_basis(2,order_))
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

  void draw() {
    // bind the buffer for the indices we want to draw
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
  }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_index> indices_;
  gl_index buffer_;
};

// forward declarations
//class TopologyBase;
//template<typename type> class Topology;

struct MeshFacet;

class Tessellation {

public:
  Tessellation( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  void build( const TopologyBase& topology );

  void draw();
  void draw_triangles( ShaderProgram& );
  void draw_edges( ShaderProgram& );

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
};

} // graphics

} // avro

#endif
