// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_TRANSPOSE_H
#define SANS_TRANSPOSE_H

#include "tools/SANSnumerics.h" //Real
#include "TransposeTraits.h"

namespace numpack 
{
// This file defines the "Transpose" function. It is needed at this
// high level due to the coupling between the different types of matrices
// that it creates. Having the Transpose functions where the transpose
// classes are defined creates touchy and complicated include patterns

//=============================================================================
//Allows the transpose to be used even when the data type is a scalar
inline int&
Transpose(int& Matrix) { return Matrix; }

inline const int&
Transpose(const int& Matrix) { return Matrix; }

inline Real&
Transpose(Real& Matrix) { return Matrix; }

inline const Real&
Transpose(const Real& Matrix) { return Matrix; }

template <int N>
inline SurrealS<N>&
Transpose(SurrealS<N>& Matrix) { return Matrix; }

template <int N>
inline const SurrealS<N>&
Transpose(const SurrealS<N>& Matrix) { return Matrix; }

//=============================================================================
//Operator to generate a datatype to represent the transpose of a matrix
template< int M, int N, class T >
inline DLA::MatrixSTranspose<M,N,T>
Transpose( const DLA::MatrixS<M,N,T>& Matrix )
{
  return DLA::MatrixSTranspose<M,N,T>( const_cast<DLA::MatrixS<M,N,T>&>( Matrix ) );
}


//=============================================================================
//Operator to generate a datatype to represent the transpose of a matrix
template< class T >
inline DLA::MatrixDTranspose< T >
Transpose(const DLA::MatrixDView<T>& Matrix)
{
  return DLA::MatrixDTranspose< T >( Matrix );
}

//Simply transpose the dimensions for a DenseNonZeroPattern
template< class T >
inline DLA::DenseNonZeroPattern< typename TransposeTraits<T>::type >
Transpose(DLA::DenseNonZeroPattern<T>& nz)
{
  return DLA::DenseNonZeroPattern< typename TransposeTraits<T>::type >( nz.n(), nz.m() );
}


//=============================================================================
//This transpose needs to be in the DLA namespace so MatrixD_Transpose can use it
template< class TM >
inline SLA::SparseMatrix_CRS_Transpose< TM >
Transpose(SLA::SparseMatrix_CRS<TM>& Matrix)
{
  return SLA::SparseMatrix_CRS_Transpose<TM>( Matrix );
}

//=============================================================================
//This transpose needs to be in the DLA namespace so MatrixD_Transpose can use it
template< class TM >
inline SLA::SparseNonZeroPattern_Transpose< TM >
Transpose( SLA::SparseNonZeroPattern<TM>& nz )
{
  return SLA::SparseNonZeroPattern_Transpose< TM >( nz );
}

//=============================================================================
//Operator to generate a datatype to represent the transpose of a matrix
template< class M00, class M01, class M10, class M11 >
inline BLA::MatrixBlock_2x2_Transpose< M00, M01, M10, M11 >
Transpose(const BLA::MatrixBlock_2x2<M00, M01, M10, M11>& Matrix)
{
  return BLA::MatrixBlock_2x2_Transpose< M00, M01, M10, M11 >( Matrix );
}

} //namespace numpack 

#endif //SANS_TRANSPOSE_H
