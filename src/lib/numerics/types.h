#ifndef URSA_LIB_NUMERICS_TYPES_H_
#define URSA_LIB_NUMERICS_TYPES_H_

#include "common/types.h"

#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/DynamicSize/VectorD.h"

namespace ursa
{

template<typename T>
class Coord : public numpack::DLA::VectorD<T>
{
protected:
  Coord( coord_t dim ) :
    numpack::DLA::VectorD<T>(dim)
  {}

  Coord( real* x , coord_t dim ) :
    numpack::DLA::VectorD<T>(dim,x)
  {}

  coord_t dim() const { return numpack::DLA::VectorD<T>::m(); }
};

class ParaCoord : public Coord<real> {};
class PhysCoord : public Coord<real> {};

class HomCoord : public Coord<real>
{
public:
  HomCoord( coord_t dim ) :
    Coord<real>( dim+1 )
  {}

  HomCoord( real* x , coord_t dim ) :
    Coord<real>( dim+1 )
  {
    for (coord_t d=0;d<dim;d++)
      this->operator()(d) = x[d];
    this->operator()(dim) = 1.0;
  }

  coord_t dim() const { return numpack::DLA::VectorD<real>::m()-1; }
};

template<typename T>
class Gradient : public numpack::DLA::VectorD<T>
{
public:
  Gradient( coord_t dim ) :
    numpack::DLA::VectorD<T>(dim)
  {}

  Gradient( real* x , coord_t dim ) :
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

} // ursa

#endif
