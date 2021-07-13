#include "element/simplex.h"

#include "graphics/new/primitives.h"
#include "graphics/new/vao.h"

#include "mesh/points.h"

namespace avro
{

namespace graphics
{

PointPrimitive::PointPrimitive( const Points& points ) :
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

PointPrimitive::~PointPrimitive() {
  glDeleteBuffers( 1 , &buffer_ );
}

void
PointPrimitive::write() {
  GL_CALL( glGenBuffers( 1 , &buffer_ ) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, buffer_ ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * coordinates_.size() , coordinates_.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
}

void
PointPrimitive::draw( ShaderProgram& program ) {
  program.use();
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , buffer_) );
  GL_CALL( glDrawArrays( GL_POINTS , 0 , coordinates_.size()/3 ) );
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, 0) );
}

void
PointPrimitive::print() const {
  for (index_t k = 0; k < nb(); k++) {
    printf("point %lu = ( ",k);
    for (coord_t d = 0; d < 3; d++)
      printf("%g ",coordinates_[3*k+d]);
    printf(")\n");
  }
}


EdgePrimitive::EdgePrimitive( coord_t order ) :
  order_(order),
  nb_basis_(order+1)
{}

void
EdgePrimitive::add( const index_t* v , index_t nv ) {
  avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
  for (index_t j = 0; j < nv; j++)
    indices_.push_back(v[j]);
}

void
EdgePrimitive::write() {
  // bind the triangles
  GL_CALL( glGenBuffers( 1 , &buffer_ ) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * indices_.size() , indices_.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
}

void
EdgePrimitive::draw() {
  // bind the buffer for the indices we want to draw
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
  GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
  GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
}

void
EdgePrimitive::print() const {
  for (index_t k = 0; k < nb(); k++) {
    printf("edge %lu = ( ",k);
    for (index_t j = 0; j < nb_basis_; j++)
      printf("%d ",int(indices_[nb_basis_*k+j]));
    printf(")\n");
  }
}

TrianglePrimitive::TrianglePrimitive( coord_t order ) :
  order_(order),
  nb_basis_(nb_simplex_basis(2,order_)),
  visible_(true)
{}

TrianglePrimitive::~TrianglePrimitive() {
  glDeleteBuffers( 1 , &buffer_ );
}

void
TrianglePrimitive::add( const index_t* v , index_t nv ) {
  avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
  for (index_t j = 0; j < nv; j++)
    indices_.push_back(v[j]);
}

void
TrianglePrimitive::write() {
  // bind the triangles
  GL_CALL( glGenBuffers( 1 , &buffer_ ) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(gl_index) * indices_.size() , indices_.data() , GL_STATIC_DRAW) );
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) );
}

void
TrianglePrimitive::draw() {

    // bind the buffer for the indices we want to draw
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
}

void
TrianglePrimitive::print() const {
  for (index_t k = 0; k < nb(); k++) {
    printf("triangle %lu = ( ",k);
    for (index_t j = 0; j < nb_basis_; j++)
      printf("%d ",int(indices_[nb_basis_*k+j]));
    printf(")\n");
  }
}

FieldData::FieldData( coord_t order ) :
  order_(order),
  nb_basis_(nb_simplex_basis(2,order_))
{}

void
FieldData::add( real_t* f , index_t ndof ) {
  avro_assert( ndof == nb_basis_ );
  for (index_t j = 0; j < ndof; j++)
    data_.push_back(f[j]);
}

FieldPrimitive::FieldPrimitive(bool nowrite) {
  // initialize the buffer
  if (!nowrite) {
    GL_CALL( glGenBuffers( 1 , &buffer_ ) );
    GL_CALL( glGenTextures( 1 , &texture_) );
  }
}

FieldPrimitive::~FieldPrimitive() {
  glDeleteBuffers( 1 , &buffer_ );
  glDeleteTextures( 1 , &texture_ );
}


void
FieldPrimitive::add( const std::string& name , index_t rank , std::shared_ptr<FieldData> data ) {
  data_.insert( { field_name(name,rank) ,data } );
}

void
FieldPrimitive::write() {

  FieldData* field = data_[active_name()].get();

  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , buffer_ ) );
  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(gl_float) * field->data().size() , field->data().data() , GL_STATIC_DRAW) );

  GL_CALL( glActiveTexture( GL_TEXTURE0 + 0 ) ); // fields are always in texture 0
  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture_) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , buffer_ ) );
}

void
FieldPrimitive::activate( ShaderProgram& shader ) {

  // bind the desired solution texture to texture unit 0
  // solution fields are alway in texture unit 0 in the shaders
  glActiveTexture(GL_TEXTURE0 + 0);

  GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture_) );
  GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , buffer_ ) );

  GLint solution_location = glGetUniformLocation(shader.handle() , "solution");
  glUniform1i(solution_location, 0); // first sampler in fragment shader
}


} // graphics

} // avro
