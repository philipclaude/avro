#ifndef LUNA_LIB_LINEAR_ALGEBRA_H_
#define LUNA_LIB_LINEAR_ALGEBRA_H_

#include "numerics/matrix.h"

#include <numpack/dense/dynamic/MatrixD_Det.h>
#include <numpack/dense/dynamic/MatrixD_Diag.h>
#include <numpack/dense/dynamic/Eigen.h>
#include <numpack/Transpose.h>

#include <math.h>

namespace luna
{

namespace numerics
{

template<typename type>
int
kernel( const MatrixD<type>& A , MatrixD<type>& K );

template<typename type>
type
determinant( const MatrixD<type>& A )
{
  return numpack::DLA::Det(A);
}

template<typename type>
type
determinant( const SymMatrixD<type>& A )
{
  return numpack::DLA::Det(A);
}

template< class T >
inline numpack::DLA::MatrixSymD<T>
log(const numpack::DLA::MatrixSymD<T>& A)
{
  // compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  // compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = std::log(LE.L[i]);

  // return the symmetric matrix
  return LE;
}

template< class T >
inline numpack::DLA::MatrixSymD<T>
exp(const numpack::DLA::MatrixSymD<T>& A)
{
  // compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  // compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = std::exp(LE.L[i]);

  // return the symmetric matrix
  return LE;
}

} // numerics

} // luna

#endif
