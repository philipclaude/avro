//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"

#include "common/types.h"

#include "numerics/vector.h"

namespace avro
{

namespace numerics
{

template<typename T>
Vector<T>
Vector<T>::operator+( const Vector<T>& Y ) const
{
  avro_assert( nb() == Y.nb() );
  Vector<T> Z( nb() );
  for (index_t k=0;k<nb();k++)
    Z[k] = data_[k] +Y[k];
  return Z;
}

template<typename T>
Vector<T>
Vector<T>::operator*( const T& a ) const
{
  Vector<T> Y( nb() );
  for (index_t k=0;k<nb();k++)
    Y[k] = a*data_[k];
  return Y;
}

template class Vector<real_t>;

} // numerics

} // avro
