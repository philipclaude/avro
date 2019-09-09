// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_OR_T_H
#define MATRIXS_OR_T_H

#include "numpack/DenseLinAlg/StaticSize/MatrixS_Type.h"

namespace numpack 
{
namespace DLA
{

template<int m, int n, class T>
struct MatrixS_or_T
{
  typedef DLA::MatrixS<m, n, T> type;
};

template<class T>
struct MatrixS_or_T<1,1,T>
{
  typedef T type;
};

template<int m, class T>
struct VectorS_or_T
{
  typedef DLA::VectorS<m, T> type;
};

template<class T>
struct VectorS_or_T<1,T>
{
  typedef T type;
};

} // namespace DLA
} // namespace numpack 

#endif // MATRIXS_OR_T_H