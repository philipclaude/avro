//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_SRC_LIB_MASTER_REFERENCE_H_
#define avro_SRC_LIB_MASTER_REFERENCE_H_

#include "avro_types.h"

#include "common/error.h"

#include "element/basis.h"

#include "numerics/functions.h"

#include <memory>
#include <vector>

namespace avro
{

class Simplex;
class Polytope;
template<typename type> class QuadratureStore;

class ReferenceElementBase
{

protected:
  ReferenceElementBase( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

public:

  coord_t order() const { return order_; }
  coord_t number() const { return number_; }

  real_t unit_volume() const { return unit_volume_; }
  real_t orthogonal_volume() const { return orthogonal_volume_; }

  index_t nb_basis() const { return nb_basis_; }
  index_t nb_barycentric() const { return number_ + 1; }

  index_t nb_interior() const { return interior_.size(); }
  index_t interior( index_t k ) const { return interior_[k]; }

  const real_t* get_reference_coordinate( index_t k ) const;
  index_t find_reference_index( const real_t* x ) const;

protected:
  coord_t number_;
  coord_t order_;
  index_t nb_basis_;

  // volume of reference element
  real_t unit_volume_;
  real_t orthogonal_volume_;

  std::vector<real_t>  nodes_;    // high-order node coordinates in canonical order
  std::vector<index_t> interior_; // indices of interior nodes
};

template<typename type>
class ReferenceElementType : public ReferenceElementBase {
public:
  ReferenceElementType( coord_t number , coord_t order ) :
    ReferenceElementBase(number,order)
  {}

  void set_basis( BasisFunctionCategory category );
  const Basis<type>& basis() const { avro_assert(basis_ != nullptr); return *basis_.get(); }
  Basis<type>& basis() { avro_assert(basis_ != nullptr); return *basis_.get(); }
  const QuadratureStore<type>& quadrature() const { avro_assert(quadrature_ != nullptr); return *quadrature_; }

protected:
  std::shared_ptr<Basis<type>> basis_;
  QuadratureStore<type>* quadrature_;
};


template<typename type> class ReferenceElement;

template<>
class ReferenceElement<Simplex> : public ReferenceElementType<Simplex> {
public:
  ReferenceElement( coord_t number , coord_t order ) :
    ReferenceElementType(number,order) {
    build();
  }

private:
  void build();
};

} // avro

#endif
