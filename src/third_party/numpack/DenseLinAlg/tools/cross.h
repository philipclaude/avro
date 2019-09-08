// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DENSELINALG_CROSS_H
#define DENSELINALG_CROSS_H

#include "tools/SANSnumerics.h"

#include "numpack/DenseLinAlg/StaticSize/MatrixS_Type.h"

#include "PromoteSurreal.h"

namespace numpack 
{

template<class T1, class T2>
inline typename promote_Surreal< T1, T2 >::type
cross(const DLA::MatrixS<2,1,T1>& a, const DLA::MatrixS<2,1,T2>& b)
{
  typename promote_Surreal< T1, T2 >::type c;

  c = a(0,0)*b(1,0) - a(1,0)*b(0,0);

  return c;
}

template<class T1, class T2>
inline DLA::MatrixS<3,1, typename promote_Surreal< T1, T2 >::type >
cross(const DLA::MatrixS<3,1,T1>& a, const DLA::MatrixS<3,1,T2>& b)
{
  DLA::MatrixS<3,1, typename promote_Surreal< T1, T2 >::type > c;

  c(0,0) = a(1,0)*b(2,0) - a(2,0)*b(1,0);
  c(1,0) = a(2,0)*b(0,0) - a(0,0)*b(2,0);
  c(2,0) = a(0,0)*b(1,0) - a(1,0)*b(0,0);

  return c;
}

}

#endif //DENSELINALG_CROSS_H
