// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DLA_INDEX_H
#define DLA_INDEX_H

#include "numpack/types/SurrealS_Type.h"
#include "numpack/DenseLinAlg/StaticSize/VectorS.h"

namespace numpack 
{
namespace DLA
{
// index is a function that gives access to an entry of a matrix/vector by index,
// or returns the object itself if it's not a matrix/vector (i.e. it's a scalar)

// scalar
inline       int& index(       int& v, const int ) { return v; }
inline const int& index( const int& v, const int ) { return v; }

inline       Real& index(       Real& v, const int ) { return v; }
inline const Real& index( const Real& v, const int ) { return v; }

inline       Real& index(       Real& m, const int, const int ) { return m; }
inline const Real& index( const Real& m, const int, const int ) { return m; }

template<int N, class T>       SurrealS<N,T>& index(       SurrealS<N,T>& v, const int ) { return v; }
template<int N, class T> const SurrealS<N,T>& index( const SurrealS<N,T>& v, const int ) { return v; }

template<int N, class T>       SurrealS<N,T>& index(       SurrealS<N,T>& m, const int, const int ) { return m; }
template<int N, class T> const SurrealS<N,T>& index( const SurrealS<N,T>& m, const int, const int ) { return m; }

// vector
template<int M, class T>       T& index(       VectorS<M,T>& v, const int n ) { return v[n]; }
template<int M, class T> const T& index( const VectorS<M,T>& v, const int n ) { return v[n]; }

template<int M, class T>       T& index(       VectorS<M,T>& v, const int i, const int ) { return v[i]; }
template<int M, class T> const T& index( const VectorS<M,T>& v, const int i, const int ) { return v[i]; }

// matrix
template<int M, int N, class T>       T& index(       MatrixS<M,N,T>& m, const int i, const int j ) { return m(i,j); }
template<int M, int N, class T> const T& index( const MatrixS<M,N,T>& m, const int i, const int j ) { return m(i,j); }


// A vector of vector effectively reverses the index relative to a matrix, hence the transposed indexing below
template<int M, int N, class T>       T& index(       VectorS<N,VectorS<M,T>>& v, const int i, const int j ) { return v[j][i]; }
template<int M, int N, class T> const T& index( const VectorS<N,VectorS<M,T>>& v, const int i, const int j ) { return v[j][i]; }

}
}

#endif //DLA_INDEX_H
