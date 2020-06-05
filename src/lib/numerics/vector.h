//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_VECTOR_H_
#define avro_LIB_NUMERICS_VECTOR_H_

#include "common/types.h"

#include <vector>

namespace avro
{

namespace numerics
{

template<typename T>
class Vector
{
public:
  Vector( const std::vector<T>& x ) :
    data_(x)
  {}

  Vector( index_t _nb ) :
    data_(_nb,T(0))
  {}

  T norm() const;
  T dot( const Vector<T>& Y ) const;

  index_t nb() const { return data_.size(); }

  T& operator[](const index_t k)
    { return data_[k]; }
  const T& operator[](const index_t k) const
    { return data_[k]; }

  Vector<T> operator*(const T& a) const;
  Vector<T>& operator*=(const T& a);

  Vector<T> operator+(const Vector<T>& Y) const;
  Vector<T>& operator+=(const Vector<T>& Y);

  Vector<T> operator-(const Vector<T>& Y) const
    { return operator+(Y*-1); }


private:
  std::vector<T> data_;

};

} // numerics

} // avro

#endif
