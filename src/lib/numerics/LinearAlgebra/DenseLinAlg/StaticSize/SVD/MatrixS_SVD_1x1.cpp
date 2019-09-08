// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MatrixS_SVD.h"
#include "../VectorS.h"
#include "../MatrixSymS.h"

#include "tools/minmax.h"
#include "tools/SANSnumerics.h"

#include "types/SurrealS.h"

namespace SANS
{
namespace DLA
{

template< int M, int N, class T >
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT )
{
  //Trivial
  U(0,0) = 1;
  S[0] = A(0,0);
  VT(0,0) = 1;
}

#define INSTANTIATE_SVD(T) \
template void SVD<1,1,T>(const MatrixS<1,1,T>& A, MatrixS<1,1,T>& U, VectorS<1,T>& S, MatrixS<1,1,T>& VT);

INSTANTIATE_SVD(Real)
INSTANTIATE_SVD(SurrealS<1>)

}
}
