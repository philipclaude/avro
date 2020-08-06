// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymD.h"
#include "../VectorD.h"

#include "MatrixSymD_Eigen_1x1.h"
#include "MatrixSymD_Eigen_2x2.h"
#include "MatrixSymD_Eigen_3x3.h"
#include "MatrixSymD_Eigen_Jacobi.h"

#include "tinymat/types/SurrealD.h"
#include "tinymat/types/SurrealS.h"

#include "tools/SANSnumerics.h"

#include <cassert>

namespace tinymat
{
namespace DLA
{

template<class T>
void
EigenValues( const MatrixSymD<T>& A, VectorD<T>& L )
{
  MatrixD<T> E(A.m(),A.n());
  if (A.m()==1)
  {
    EigenSystem_1x1(A,L,E);
  }
  else if (A.m()==2)
  {
    EigenSystem_2x2(A,L,E);
  }
  else if (A.m()==3)
  {
    EigenSystem_Jacobi(A,L,E);
  }
  else
    EigenSystem_Jacobi(A,L,E);
}

template<class T>
void
EigenVectors( const MatrixSymD<T>& A, MatrixD<T>& E )
{
  VectorD<T> L(A.m());
  if (A.m()==1)
  {
    EigenSystem_1x1(A,L,E);
  }
  else if (A.m()==2)
  {
    EigenSystem_2x2(A,L,E);
  }
  else if (A.m()==3)
  {
    EigenSystem_Jacobi(A,L,E);
  }
  else
    EigenSystem_Jacobi(A,L,E);
}

template<class T>
void
EigenSystem( const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E )
{
  if (A.m()==1)
  {
    EigenSystem_1x1(A,L,E);
  }
  else if (A.m()==2)
  {
    EigenSystem_2x2(A,L,E);
  }
  else if (A.m()==3)
  {
    EigenSystem_Jacobi(A,L,E);
  }
  else
    EigenSystem_Jacobi(A,L,E);
}

#define INSTANTIATE_EIGEN(T) \
template void EigenSystem<T>(const MatrixSymD<T>& A, VectorD<T>& L, MatrixD<T>& E ); \
template void EigenValues<T>(const MatrixSymD<T>& A, VectorD<T>& L ); \
template void EigenVectors<T>(const MatrixSymD<T>& A, MatrixD<T>& L );

INSTANTIATE_EIGEN(Real)
INSTANTIATE_EIGEN(SurrealS<1>)
INSTANTIATE_EIGEN(SurrealS<3>)
INSTANTIATE_EIGEN(SurrealS<6>)
INSTANTIATE_EIGEN(SurrealS<10>)


}
}
