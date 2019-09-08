// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_SVD_H
#define MATRIXS_SVD_H

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "tools/minmax.h"
#include "tools/SANSException.h"

//SVD for DenseLinAlg.

namespace SANS
{
namespace DLA
{

// Singular value decomposition such that
// A = U * S * VT
template< int M, int N, class T >
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT );

template< int M, class T >
void
SVD( const MatrixSymS<M,T>& A, MatrixS<M,M,T>& U, VectorS<M,T>& S, MatrixS<M,M,T>& VT )
{
  MatrixS<M,M,T> Amat = A;
  SVD(Amat, U, S, VT);
}

} //namespace SANS
} //namespace DLA

#endif //MATRIXS_SVD_H
