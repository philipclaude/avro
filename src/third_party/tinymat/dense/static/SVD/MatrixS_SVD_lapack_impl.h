// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXS_SVD_LAPACK_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "tools/SANSException.h"

#ifdef DLA_LAPACK

#include "../MatrixS_SVD.h"
#include "tinymat/dense/static/MatrixS.h"
#include "tinymat/dense/static/VectorS.h"
#include "tinymat/dense/tools/dense_LAPACK.h"

//Singular Value Decomposition. The lapack routines for single/double are defined
//as part of the explicit instantiation in the .cpp files

namespace tinymat 
{
namespace DLA
{

// A is the matrix that is decomposed
// U and VT are orthogonal matricies, and S is a diagonal matrix
// Note that VT is Transpose(V), hence
// A = U * S * VT
template<int M, int N, class T>
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT )
{
  int m = M;
  int n = N;
  char jobu = 'A'; //Compute all of U matrix
  char jobv = 'A'; //Compute all of V matrix
  int INFO;
  T work[MAX(3*MIN(M,N) + MAX(M,N), 5*MIN(M,N))]; // Work array
  int lwork = sizeof(work)/sizeof(T); // Work array size
  MatrixS<M,N,T> Acp(A); // Copy of the matrix

  //Note that VT and U are swapped (as well as m,n) for the transpose of C vs Fortran ordering
  LAPACK_GESVD(&jobu, &jobv, &n, &m, &Acp(0,0), &n, &S[0], &VT(0,0), &n, &U(0,0), &m, work, &lwork, &INFO);

  SANS_ASSERT_MSG( INFO == 0, "INFO == %d", INFO );
}

#define INSTANTIATE_SVD(M,N,T) \
template void SVD<M,N,T>(const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT);

} //namespace tinymat 
} //namespace DLA

#endif
