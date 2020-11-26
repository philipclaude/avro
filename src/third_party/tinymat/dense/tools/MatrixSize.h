// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SRC_LINEARALGEBRA_dense_TOOLS_MATRIXSIZE_H_
#define SRC_LINEARALGEBRA_dense_TOOLS_MATRIXSIZE_H_

#include "tools/SANSnumerics.h"
#include "tinymat/dense/static/MatrixS_Type.h"

// forward declaration
template<int N, class T>
class SurrealS;

namespace tinymat 
{
namespace DLA
{

template<class T>
struct MatrixSize;

// Specializations: Gives the size of a MatrixS or Real
template<int M_, int N_, class T>
struct MatrixSize<MatrixS<M_,N_,T> >
{
  static const int M = M_;
  static const int N = N_;
};

template<>
struct MatrixSize<Real>
{
  static const int M = 1;
  static const int N = 1;
};

template<int N_, class T>
struct MatrixSize<SurrealS<N_,T>>
{
  static const int M = 1;
  static const int N = 1;
};

}
}

#endif /* SRC_LINEARALGEBRA_dense_TOOLS_MATRIXSIZE_H_ */
