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
  nb_basis_(order+1),
  visible_(true)
{}

void
EdgePrimitive::add( const index_t* v , index_t nv ) {
  avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
  for (index_t j = 0; j < nv; j++)
    indices_.push_back(v[j]);
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

void
TrianglePrimitive::add( const index_t* v , index_t nv ) {
  avro_assert_msg( nv == nb_basis_ , "nv = %lu, nb_basis = %lu" , nv , nb_basis_ );
  for (index_t j = 0; j < nv; j++)
    indices_.push_back(v[j]);
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

void
FieldPrimitive::add( const std::string& name , index_t rank , std::shared_ptr<FieldData> data ) {
  data_.insert( { field_name(name,rank) ,data } );
}

} // graphics

} // avro
