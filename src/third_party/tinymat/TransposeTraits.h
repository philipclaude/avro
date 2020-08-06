// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_TRANSPOSETRAITS_H
#define SANS_TRANSPOSETRAITS_H

#include "tinymat/dense/dynamic/MatrixD_Type.h"
#include "tinymat/dense/static/MatrixS_Type.h"
#include "tinymat/types/SurrealS_Type.h"

namespace tinymat
{

//=============================================================================
// Defines transpose view types, such as MatrixSTranspose and SparseMatrix_CRS_Transpose
template<class T>
struct TransposeViewTraits;

// Defines actual transpose of a matrix, i.e. TransposeTraits<MatrixS<M,N>>::type is MatrixS<N,M>
template<class T>
struct TransposeTraits;

//=============================================================================
// Scalar transposes
template<>
struct TransposeViewTraits< int >
{
  typedef int& type;
};

template<>
struct TransposeTraits< int >
{
  typedef int type;
};

template<>
struct TransposeViewTraits< float >
{
  typedef float& type;
};

template<>
struct TransposeTraits< float >
{
  typedef float type;
};

template<>
struct TransposeViewTraits< double >
{
  typedef double& type;
};

template<>
struct TransposeTraits< double >
{
  typedef double type;
};

template<int N>
struct TransposeViewTraits< SurrealS<N> >
{
  typedef SurrealS<N>& type;
};

template<int N>
struct TransposeTraits< SurrealS<N> >
{
  typedef SurrealS<N> type;
};


//=============================================================================
// Statically sized matrix
template<int M, int N, class T>
struct TransposeViewTraits< DLA::MatrixS<M,N,T> >
{
  typedef DLA::MatrixSTranspose<M,N,T> type;
};

template<int M, int N, class T>
struct TransposeTraits< DLA::MatrixS<M,N,T> >
{
  typedef DLA::MatrixS<N,M, typename TransposeTraits<T>::type > type;
};

template<int M, class T>
struct TransposeViewTraits< DLA::VectorS<M,T> >
{
  typedef DLA::MatrixSTranspose<M,1,T> type;
};

template<int M, class T>
struct TransposeTraits< DLA::VectorS<M,T> >
{
  typedef DLA::MatrixS<1,M, typename TransposeTraits<T>::type > type;
};

//=============================================================================
// Dynamically sized matrix
template<class T>
struct TransposeViewTraits< DLA::MatrixD<T> >
{
  typedef DLA::MatrixDTranspose< T > type;
};

template<class T>
struct TransposeViewTraits< DLA::MatrixDView<T> >
{
  typedef DLA::MatrixDTranspose< T > type;
};

template<class T>
struct TransposeTraits< DLA::MatrixD<T> >
{
  typedef DLA::MatrixD< typename TransposeTraits<T>::type > type;
};


//=============================================================================
// Dense dummy non-zero pattern transpose
template<class T>
struct TransposeViewTraits< DLA::DenseNonZeroPattern<T> >
{
  typedef DLA::DenseNonZeroPattern< typename TransposeTraits<T>::type > type;
};

template<class T>
struct TransposeTraits< DLA::DenseNonZeroPattern<T> >
{
  typedef DLA::DenseNonZeroPattern< typename TransposeTraits<T>::type > type;
};

} // namespace tinymat


#endif // SANS_TRANSPOSETRAITS_H
