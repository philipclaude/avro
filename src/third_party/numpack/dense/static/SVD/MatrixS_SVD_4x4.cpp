// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MatrixS_SVD.h"
#include "../VectorS.h"
#include "../MatrixSymS.h"

#include "tools/SANSnumerics.h"

#include "numpack/types/SurrealS.h"

namespace numpack 
{
namespace DLA
{

#ifndef DLA_LAPACK

template< int M, int N, class T >
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT )
{
  SANS_DEVELOPER_EXCEPTION("not implemented -- turn LAPACK on");
}


#define INSTANTIATE_SVD(T) \
template void SVD<4,4,T>(const MatrixS<4,4,T>& A, MatrixS<4,4,T>& U, VectorS<4,T>& S, MatrixS<4,4,T>& VT);

INSTANTIATE_SVD(Real)

#endif

}
}
