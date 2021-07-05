//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_SRC_LIB_MASTER_REFERENCE_H_
#define avro_SRC_LIB_MASTER_REFERENCE_H_

#include "avro_types.h"

#include "numerics/functions.h"

#include <vector>

namespace avro
{

class Simplex;
class Polytope;
class Hypercube;

template<typename type, int N, int P> struct LagrangeNodes;

template<int P>
struct LagrangeNodes<Simplex,1,P> {
  static const std::vector<real_t> coord_s_,coord_t_;
};

template<int P>
struct LagrangeNodes<Simplex,2,P> {
  static const std::vector<real_t> coord_s_,coord_t_;
};

template<int P>
struct LagrangeNodes<Simplex,3,P> {
  static const std::vector<real_t> coord_s_,coord_t_,coord_u_;
};

template<int P>
struct LagrangeNodes<Simplex,4,P> {
  static const std::vector<real_t> coord_s_,coord_t_,coord_u_,coord_v_;
};

class ReferenceElementBase
{
public:
  ReferenceElementBase( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

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

template<typename type> class ReferenceElement;

template<typename type>
class ReferenceElement {

public:
  ReferenceElement( coord_t number , coord_t order );

  const real_t* get_reference_coordinate( index_t k ) const;

  real_t vunit() const { return unit_volume_; }
  index_t nb_basis() const;

  index_t nb_interior() const { return interior_.size(); }
  index_t interior( index_t k ) const { return interior_[k]; }

  const index_t* get_lattice_coordinate( index_t k ) const;
  int find_index( const index_t* x ) const;
  int find_index( const real_t* x ) const;

private:
  coord_t number_;
  coord_t order_;
  std::vector<real_t> nodes_;
  std::vector<index_t> lattice_;
  std::vector<index_t> interior_;

  real_t unit_volume_;
  real_t orth_volume_;
};

/*
template<>
class ReferenceElement<Simplex> : public ReferenceElementBase
{
public:
  ReferenceElement( coord_t number , coord_t order ) :
    ReferenceElementBase(number,order)
  {
    //precalculate();
  }

  const real_t* get_reference_coordinate( index_t k ) const;
  //const index_t* get_lattice_coordinate( index_t k ) const;

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
  //std::vector<index_t> interior_;

  std::vector<real_t> nodes_;

};

template<>
class ReferenceElement<Polytope> : public ReferenceElementBase
{
public:
  ReferenceElement( coord_t number , coord_t order ) :
    ReferenceElementBase(number,order)
  {}

};
*/

} // avro

#endif
