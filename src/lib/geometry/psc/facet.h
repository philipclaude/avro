//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_PSC_FACET_H_
#define avro_LIB_GEOMETRY_PSC_FACET_H_

#include "geometry/psc/object.h"
#include "numerics/mat.h"
#include "numerics/vec.h"

namespace avro
{

namespace PSC
{

class Facet : public Object
{
public:
  Facet( Body* body , std::vector<std::shared_ptr<Entity>>& facets );
  void build_basis();

  ~Facet() {}

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const
    { inverse(x,u); }
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const;

  const matd<real_t>& basis() const { return V_; }
  int dimension() const { return dimension_; }

  void set_basis_by_name();
  const vecd<real_t>& x0() const { return x0_; }

private:
  matd<real_t> V_;
  matd<real_t> B_;
  vecd<real_t> x0_;

  int dimension_;
};

} // PSC

} // avro

#endif
