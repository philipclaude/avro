//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MASTER_QUADRATURE_H_
#define avro_LIB_MASTER_QUADRATURE_H_

#include "common/types.h"

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

} // avro

#endif
