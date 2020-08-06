//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_TYPES_H_
#define avro_LIB_NUMERICS_TYPES_H_

#include "common/types.h"

#include "tinymat/dense/dynamic/MatrixD.h"
#include "tinymat/dense/dynamic/VectorD.h"

namespace avro
{

template<typename T>
class Coord : public tinymat::DLA::VectorD<T>
{
protected:
  Coord( coord_t dim ) :
    tinymat::DLA::VectorD<T>(dim)
  {}

  Coord(real_t* x , coord_t dim ) :
    tinymat::DLA::VectorD<T>(dim,x)
  {}

  coord_t dim() const { return tinymat::DLA::VectorD<T>::m(); }
};

class ParaCoord : public Coord<real_t> {};
class PhysCoord : public Coord<real_t> {};

class HomCoord : public Coord<real_t>
{
public:
  HomCoord( coord_t dim ) :
    Coord<real_t>( dim+1 )
  {}

  HomCoord(real_t* x , coord_t dim ) :
    Coord<real_t>( dim+1 )
  {
    for (coord_t d=0;d<dim;d++)
      this->operator()(d) = x[d];
    this->operator()(dim) = 1.0;
  }

  coord_t dim() const { return tinymat::DLA::VectorD<real_t>::m()-1; }
};

template<typename T>
class Gradient : public tinymat::DLA::VectorD<T>
{
public:
  Gradient( coord_t dim ) :
    tinymat::DLA::VectorD<T>(dim)
  {}

  Gradient(real_t* x , coord_t dim ) :
    tinymat::DLA::VectorD<T>(dim,x)
  {}

  coord_t dim() const { return tinymat::DLA::VectorD<T>::m(); }
};

template<typename T>
class Hessian : public tinymat::DLA::MatrixD<T>
{
public:
  Hessian( coord_t dim ) :
    tinymat::DLA::MatrixD<T>(dim,dim)
  {}
};

} // avro

#endif
