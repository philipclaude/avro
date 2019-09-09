// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymD.h"
#include "../VectorD.h"

#include "../MatrixD_Det.h"
#include "../MatrixD_Trace.h"

#include "numpack/types/SurrealD.h"

#include "tools/SANSnumerics.h"

#include <cassert>

namespace numpack
{
namespace DLA
{


template<class T>
void
EigenValues(const MatrixSymD<T>& A, VectorD<T>& L )
{
  //Trivial
  L[0] = A(0,0);
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

  EigenValues( A, L );

  if (A.m()==1)
  {
    //Trivial
    E(0,0) = 1;
  }
  else
    assert(false);
}

#define INSTANTIATE_EIGEN(T) \
template void EigenValues<T>(const MatrixSymD<T>& A, VectorD<T>& L ); \
template void EigenVectors<T>(const MatrixSymD<T>& A, MatrixD<T>& E ); \
template void EigenSystem<T>(const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E );

INSTANTIATE_EIGEN(Real)
INSTANTIATE_EIGEN(SurrealD)

//INSTANTIATE_EIGEN(SurrealS<1>)
//INSTANTIATE_EIGEN(SurrealS<2>)

}
}
