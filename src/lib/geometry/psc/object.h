//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_PLC_H_
#define avro_LIB_GEOMETRY_PLC_H_

#include "avro_types.h"

#include "geometry/body.h"
#include "geometry/entity.h"

#include "numerics/mat.h"

namespace avro
{

namespace PSC
{

class Object : public Entity
{
public:
  Object( coord_t num , coord_t dim );

  virtual void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const = 0;
  virtual void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const = 0;
  virtual void evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const = 0;

  void tessellate( BodyTessellation& body_tess ) const;

  void print(bool with_children) const;
  void build_hierarchy();

protected:
  coord_t dim_;
  matd<real_t> basis_;
};

class Body : public avro::Body
{
public:
  Body( coord_t number , coord_t dim ) :
    avro::Body(number),
    dim_(dim)
  {}

  ~Body() {}

  coord_t dim() const { return dim_; }

private:
  coord_t dim_;
};

} // PSC

} // avro

#endif
