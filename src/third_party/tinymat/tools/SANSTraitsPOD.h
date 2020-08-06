// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANSTRAITSPOD_H
#define SANSTRAITSPOD_H

#include "tools/SANSnumerics.h"   // Real
#include "tools/SANSException.h"

namespace tinymat 
{

// Traits class to allow templated classes to allow Real arguments in arithmetic
// operations; examples of classes that use this are DenseMatrix<MatrixS<N,Real>>
// and VectorS<N, SurrealS<M>>

struct SANSDummyType
{
  operator int() const { SANS_ASSERT(false); return 0; }
};

template< class T >
struct POD { typedef Real type; };

template<>
struct POD<Real> { typedef SANSDummyType type; };

template<>
struct POD<int> { typedef SANSDummyType type; };

template<>
struct POD<unsigned int> { typedef SANSDummyType type; };

}  // namespace tinymat 

#endif  // SANSTRAITSPOD_H
