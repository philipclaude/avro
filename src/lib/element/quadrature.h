//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MASTER_QUADRATURE_H_
#define avro_LIB_MASTER_QUADRATURE_H_

#include "avro_types.h"

#include "element/basis.h"
#include "numerics/mat.h"

#include <map>
#include <memory>
#include <vector>

namespace avro
{

class Quadrature
{
public:
  Quadrature( coord_t dim , const int order=-1 ); // default to max quad
  void retrieve( std::vector<real_t>& x , std::vector<real_t>& w )
  {
    x = x_;
    w = w_;
  }
  virtual void define() = 0;

  coord_t dim() const { return dim_; }
  index_t nb() const { return nb_quad_; }
  int order() const { return order_; }
  int& order() { return order_; }
  bool defined() const { return defined_; }

  const real_t* x( index_t k ) const;
  real_t  w( index_t k ) const;

  real_t volume() const { return volume_; }

protected:
  virtual ~Quadrature() {}

  coord_t dim_;
  int order_;
  index_t nb_quad_;
  real_t volume_;

  bool defined_;

  std::vector<real_t> x_;
  std::vector<real_t> w_;
};

class ConicalProductQuadrature : public Quadrature
{
public:
  using Quadrature::Quadrature;

  void define();

  const coord_t order_max_ = 8;
};

class GrundmannMoellerQuadrature :   public Quadrature
{
public:
  using Quadrature::Quadrature;

  void define();

private:
  const coord_t order_max_ = 14;
  index_t rule_;
};

template<typename type>
class QuadratureStore
{
public:
  QuadratureStore( BasisFunctionCategory category );

public:
  const matd<real_t>& get_basis( index_t n , index_t p , int q=-1 ) const {
    if (q < 0) q = qmax-1;
    return basis_.at(n).at(p).at(q)[0];
  }

  const matd<real_t>& get_deriv( index_t n , index_t p , index_t q , coord_t dim ) const {
    avro_assert( dim < basis_.at(n).at(p).at(q).size()-1 );
    return basis_.at(n).at(p).at(q)[dim+1];
  }

  const Quadrature& quadrature( index_t n , int q=-1 ) const {
    if (q < 0) return *quadrature_.at(n).at(qmax-1);
    return *quadrature_.at(n).at(q);
  }

  index_t memory() const { return memory_; } // in bytes

  index_t max_quad() const { return qmax-1; }

private:
  void build();

  index_t nmax, pmax, qmax;

  const BasisFunctionCategory category_;

  // map from [shape number] -> { map from [basis function order] -> { map from [quadrature order] -> {basis,basis_s,basis_t,basis_u,basis_v} }
  std::map< coord_t , std::map<coord_t, std::map<coord_t,std::vector<matd<real_t>>>> > basis_;
  std::map< coord_t , std::map< coord_t , std::shared_ptr<Quadrature> > > quadrature_;

  index_t memory_;

public:
  /*static QuadratureStore<type>* get_lagrange() {
    static QuadratureStore<type> instance(BasisFunctionCategory_Lagrange);
    return &instance;
  }*/
};

class QuadraturePoint {
public:
  QuadraturePoint( const Quadrature& quad ) :
    quadrature_(quad),
    idx_(quad.nb())
  {}

  void set_index( index_t idx ) { idx_ = idx; }
  index_t index() const { return idx_; }

  const real_t* coordinate() const { return quadrature_.x(idx_); }
  real_t weight() const { return quadrature_.w(idx_); }

private:
  const Quadrature& quadrature_;
  index_t idx_;
};

extern QuadratureStore<Simplex> __store_simplex_legendre__;
extern QuadratureStore<Simplex> __store_simplex_lagrange__;
extern QuadratureStore<Simplex> __store_simplex_bernstein__;

} // avro

#endif
