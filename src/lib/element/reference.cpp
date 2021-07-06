//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"
#include "common/tools.h"

#include "element/lagrange_nodes.h"
#include "element/reference.h"
#include "element/quadrature.h"
#include "element/simplex.h"

#include "numerics/functions.h"
#include "numerics/geometry.h"
#include "numerics/mat.h"

#include <cmath>

namespace avro
{

void
ReferenceElement<Simplex>::build() {

  // number of barycentric coordinates
  index_t nb_coord    = nb_barycentric();
  index_t nb_vertices = number_ +1;

  // volume properties
  unit_volume_       = sqrt(number_+1)/(numerics::factorial(number_)*sqrt(pow(2.,number_)));
  orthogonal_volume_ = 1./numerics::factorial(number_);

  // compute the coordinates of the lagrange nodes
  nb_basis_ = nb_simplex_basis(number_,order_);
  if (order_ > 0) get_lagrange_nodes<Simplex>(number_,order_,nodes_);
  else {
    // yes, p = 0 Lagrange doesn't make sense, but it's okay
    for (coord_t d = 0; d < nb_coord; d++)
      nodes_.push_back( 1./nb_vertices );
  }
  avro_assert_msg( nodes_.size() == nb_coord*nb_basis_ ,
                   "|nodes| = %lu, nb_basis = %lu, number = %u", nodes_.size(),nb_basis_,number_ );

  // search for indices of interior nodes
  // interior nodes have lattice coordinates / canonical reference coordinates that are all nonzero
  // yes, this uses floating-point checks, but we shouldn't expect precision issues here
  // since this is the canonical reference elements (well spaced node coordinates)
  for (index_t k = 0; k < nb_basis_; k++) {
    bool anyzero = false;
    for (coord_t j = 0; j < nb_coord; j++) {
      if (fabs(nodes_[k*nb_coord+j]) < 1e-12 ) {
        anyzero = true;
        break;
      }
    }
    if (!anyzero) // none of the coordinates are zero, so it is an interior node
      interior_.push_back(k);
  }
  avro_assert_msg( interior_.size() == nb_simplex_basis_interior(number_,order_) ,
                   "|interior| = %lu, should be %lu", interior_.size() , nb_simplex_basis_interior(number_,order_) );
}

void
ReferenceElement<Simplex>::set_basis( BasisFunctionCategory category )
{
  basis_ = std::make_shared<Basis<Simplex>>(number_,order_,category);

  if (category == BasisFunctionCategory_Lagrange)
    quadrature_ = &__store_simplex_lagrange__;
  else if (category == BasisFunctionCategory_Legendre)
    quadrature_ = &__store_simplex_legendre__;
  else if (category == BasisFunctionCategory_Bernstein)
    quadrature_ = &__store_simplex_bernstein__;
  else if (category == BasisFunctionCategory_None) {
    // nothing to do, assume the caller knows what they're doing
  }
  else {
    printf("unknown category %d\n",category);
    avro_assert_not_reached;
  }
}

index_t
ReferenceElementBase::find_reference_index( const real_t* x ) const {
  // yes, this uses floating-point checks, but we shouldn't expect precision issues here
  // since this is the canonical reference elements (well spaced node coordinates)
  const index_t nb_coord = nb_barycentric();
  for (index_t k = 0; k < nb_basis_; k++) {
    real_t d = numerics::distance( &nodes_[k*nb_coord] , x , nb_coord );
    if (d < 1e-12) return k;
  }
  std::vector<real_t> xp(x,x+nb_coord);
  printf("number = %u, order = %u\n",number_,order_);
  print_inline(xp,"could not find");
  print_inline(nodes_);
  avro_assert_not_reached;
  return nb_basis_;
}

const real_t*
ReferenceElementBase::get_reference_coordinate( const index_t k ) const {
  const index_t nb_coord = nb_barycentric();
  avro_assert( k < nodes_.size()/nb_coord );
  return &nodes_[k*nb_coord];
}

} // avro
