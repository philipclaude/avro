// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_TRANSPOSETRAITS_H
#define SANS_TRANSPOSETRAITS_H

#include "numpack/dense/dynamic/MatrixD_Type.h"
#include "numpack/dense/static/MatrixS_Type.h"
#include "numpack/sparse/sparse_Type.h"
#include "numpack/types/SurrealS_Type.h"

namespace numpack 
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


//=============================================================================
// Sparse Matrix transpose traits
template<class TM>
struct TransposeViewTraits< SLA::SparseMatrix_CRS<TM> >
{
  typedef SLA::SparseMatrix_CRS_Transpose< TM > type;
};

template<class TM>
struct TransposeTraits< SLA::SparseMatrix_CRS<TM> >
{
  typedef SLA::SparseMatrix_CRS< typename TransposeTraits<TM>::type > type;
};

//=============================================================================
// Sparse Matrix non-zero pattern transpose
template<class TM>
struct TransposeViewTraits< SLA::SparseNonZeroPattern<TM> >
{
  typedef SLA::SparseNonZeroPattern_Transpose< TM > type;
};

template<class TM>
struct TransposeTraits< SLA::SparseNonZeroPattern<TM> >
{
  typedef SLA::SparseNonZeroPattern< typename TransposeTraits<TM>::type > type;
};


//=============================================================================
// Block matrix transposed view
namespace BLA
{
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
class MatrixBlock_2x2;

template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
class MatrixBlock_2x2_Transpose;
}

template<class M00, class M01, class M10, class M11>
struct TransposeViewTraits< BLA::MatrixBlock_2x2<M00, M01, M10, M11> >
{
  typedef BLA::MatrixBlock_2x2_Transpose< M00, M01, M10, M11 > type;
};

} // namespace numpack 


#endif // SANS_TRANSPOSETRAITS_H
