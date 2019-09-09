// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixSymD_Eigen_impl.h"

#include "numpack/types/SurrealD.h"

#include "tools/SANSnumerics.h"

namespace numpack
{
namespace DLA
{

template<class T>
void
EigenValues(const MatrixSymD<T>& A, VectorD<T>& L )
{
  MatrixD<T> E( A.m() , A.m() );
  EigenSystem( A, L, E );
}

template<class T>
void
EigenVectors(const MatrixSymD<T>& A, MatrixD<T>& E )
{
  VectorD<T> L( A.m() );
  EigenSystem( A, L, E );
}

template<class T>
void
EigenSystem(const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E )
{
  // iterate with Jacobi to find the solution
  EigenSystem_Jacobi(A, L, E);
}

#define INSTANTIATE_EIGEN(T) \
template void EigenValues<T>(const MatrixSymD<T>& A, VectorD<T>& L ); \
template void EigenVectors<T>(const MatrixSymD<T>& A, MatrixD<T>& E ); \
template void EigenSystem<T>(const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E );

INSTANTIATE_EIGEN(Real)
/*
INSTANTIATE_EIGEN(SurrealS<1>)
INSTANTIATE_EIGEN(SurrealS<10>)
INSTANTIATE_EIGEN(SurrealS<20>)
INSTANTIATE_EIGEN(SurrealS<50>)
*/
}
}
