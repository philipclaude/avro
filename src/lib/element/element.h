//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_SHAPE_SHAPE_H_
#define avro_LIB_SHAPE_SHAPE_H_

#include "common/types.h"

#include "element/basis.h"
#include "element/reference.h"

//#include "numerics/matrix.h"

#include <memory>
#include <string>
#include <vector>

namespace avro
{

class Quadrature;
template<typename type> class Basis;

class ElementBase
{
public:
  ElementBase( const std::string& name ) :
    name_(name)
  {}
  const std::string name() const { return name_; }
protected:
  std::string name_;
};

template<typename type>
class Element : public ElementBase
{
public:

  Element( coord_t number , coord_t order ) :
    ElementBase("unknown"),
    number_(number),
    order_(order),
    reference_(number_,order_),
    parameter_(false),
    basis_(nullptr)
  {}

  Element( coord_t number , coord_t order , const std::string& name ) :
    ElementBase(name),
    number_(number),
    order_(order),
    reference_(number,order_),
    parameter_(false),
    basis_(nullptr)
  {}

  void set_basis( BasisFunctionCategory category );

  coord_t number() const { return number_; }
  coord_t order() const {return order_; }

  index_t nb_quad() const { return wquad_.size(); }

  void load_quadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

  void eval( const double* x , double* phi )
  {
    avro_assert( basis_!=nullptr );
    basis_->evaluate(x,phi);
  }

  void set_parameter( bool x ) { parameter_ = x; }
  bool parameter() const { return parameter_; }

  const ReferenceElement<type>& reference() const { return reference_; }

  real_t quad_weight(index_t k) const { return wquad_[k]; }
  const real_t* quad_point(index_t k) const { return &xquad_[number_*k]; }

  const Basis<type>& basis() const { avro_assert(basis_!=nullptr); return *basis_.get(); }
  Basis<type>& basis() { avro_assert(basis_!=nullptr); return *basis_.get(); }

protected:
  coord_t number_;
  coord_t order_;

  //numerics::MatrixD<real_t> phi_;
  //std::vector< numerics::MatrixD<real_t> > dphi_;

  std::vector<real_t> xquad_;
  std::vector<real_t> wquad_;

  ReferenceElement<type> reference_;

  bool parameter_;

private:
  std::shared_ptr< Basis<type> > basis_;

};

typedef struct
{
  std::vector<index_t> indices;
  coord_t dim;
  bool sorted = true;
} ElementIndices;

bool operator< ( const ElementIndices& f , const ElementIndices& g );
bool operator== ( const ElementIndices& f , const ElementIndices& g );

} // avro

#endif
