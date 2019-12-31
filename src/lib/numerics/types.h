#ifndef avro_LIB_NUMERICS_TYPES_H_
#define avro_LIB_NUMERICS_TYPES_H_

#include "common/types.h"

#include "numpack/dense/dynamic/MatrixD.h"
#include "numpack/dense/dynamic/VectorD.h"

namespace avro
{

template<typename T>
class Coord : public numpack::DLA::VectorD<T>
{
protected:
  Coord( coord_t dim ) :
    numpack::DLA::VectorD<T>(dim)
  {}

  Coord(real_t* x , coord_t dim ) :
    numpack::DLA::VectorD<T>(dim,x)
  {}

  coord_t dim() const { return numpack::DLA::VectorD<T>::m(); }
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

  coord_t dim() const { return numpack::DLA::VectorD<real_t>::m()-1; }
};

template<typename T>
class Gradient : public numpack::DLA::VectorD<T>
{
public:
  Gradient( coord_t dim ) :
    numpack::DLA::VectorD<T>(dim)
  {}

  Gradient(real_t* x , coord_t dim ) :
    numpack::DLA::VectorD<T>(dim,x)
  {}

  coord_t dim() const { return numpack::DLA::VectorD<T>::m(); }
};

template<typename T>
class Hessian : public numpack::DLA::MatrixD<T>
{
public:
  Hessian( coord_t dim ) :
    numpack::DLA::MatrixD<T>(dim,dim)
  {}
};

} // avro

#endif
