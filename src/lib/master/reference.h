#ifndef LUNA_SRC_LIB_MASTER_REFERENCE_H_
#define LUNA_SRC_LIB_MASTER_REFERENCE_H_

#include "common/types.h"

#include "numerics/functions.h"

#include <vector>

namespace luna
{

class Simplex;
class Polytope;
class Hypercube;

class ReferenceElementBase
{
public:
  ReferenceElementBase( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  // this will create a linker error if not defined by the base
  const real_t* get_reference_coordinate( index_t k ) const;
  const index_t* get_lattice_coordinate( index_t k ) const;

  coord_t order() const { return order_; }
  coord_t number() const { return number_; }

  real_t vunit() const { return vunit_; }
  real_t vorth() const { return vorth_; }


protected:
  coord_t number_;
  coord_t order_;

  // vertex coordinates
  std::vector<real_t> xunit_; // unit element coordinates (unit edge lengths)
  std::vector<real_t> xorth_; // orthogonal corner at origin

  // all coordinates
  std::vector<real_t>  xref_;
  std::vector<index_t> lref_;

  real_t vunit_;
  real_t vorth_;
};

template<typename Shape> class ReferenceElement;

template<>
class ReferenceElement<Simplex> : public ReferenceElementBase
{
public:
  ReferenceElement( coord_t number , coord_t order ) :
    ReferenceElementBase(number,order)
  {
    precalculate();
  }

  const real_t* get_reference_coordinate( index_t k ) const;
  const index_t* get_lattice_coordinate( index_t k ) const;

  index_t nb_basis() const
  {
    index_t np = 1;
    for (coord_t d=1;d<=number_;d++)
      np *= (order_+d);
    return np/numerics::factorial(number_);
  }

  int find_index( const index_t* x ) const;
  int find_index( const real_t* x ) const;

  index_t nb_interior() const { return interior_.size(); }
  index_t interior( index_t k ) const { return interior_[k]; }

  void precalculate();

private:
  std::vector<index_t> interior_;

};

template<>
class ReferenceElement<Polytope> : public ReferenceElementBase
{
public:
  ReferenceElement( coord_t number , coord_t order ) :
    ReferenceElementBase(number,order)
  {}

};

#if 0


template<>
class ReferenceElement<Hypercube> : public ReferenceElementBase
{

};

#endif

} // luna

#endif
