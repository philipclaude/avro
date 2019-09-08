// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_ARRAY2D_H
#define MATRIXD_ARRAY2D_H

#include "tools/SANSnumerics.h"
#include "MatrixD.h"

namespace SANS
{
namespace DLA
{

//-----------------------------------------------------------------------------
//
// MatrixD_Array provides a two-dimensional array of equally sized matrices.
// The implementation uses continuous memory
//
//-----------------------------------------------------------------------------

template< class T >
class MatrixDView_Array2D
{
public:
  MatrixDView_Array2D( T v0[], const int m, const int n, const int size0, const int size1 ) :
    v_(v0), m_(m), n_(n), Msize_(m*n), size0_(size0), size1_(size1)
  {
    SANS_ASSERT( m_ >= 0 );
    SANS_ASSERT( n_ >= 0 );
  }

  MatrixDView_Array2D( const MatrixDView_Array2D& A ) :
    v_(A.v_), m_(A.m_), n_(A.n_), Msize_(A.Msize_), size0_(A.size0_), size1_(A.size1_) {}
  MatrixDView_Array2D& operator=( const MatrixDView_Array2D& ) = delete;

  int size() const { return size0_*size1_; }
  int size0() const { return size0_; }
  int size1() const { return size1_; }
  int m() const { return m_; }
  int n() const { return n_; }

  // value accessor operators
        MatrixDView<T> operator()( const int i, const int j )       { return MatrixDView<T>(v_ + (i*size1_ + j)*Msize_, m_, n_); }
  const MatrixDView<T> operator()( const int i, const int j ) const { return MatrixDView<T>(v_ + (i*size1_ + j)*Msize_, m_, n_); }

  // element wise accessors
        T& value(const int i)       { return v_[i]; }
  const T& value(const int i) const { return v_[i]; }

  T operator=( const T& val )
  {
    const int size0 = size0_;
    const int size1 = size1_;
    for (int i = 0; i < size0; ++i)
      for (int j = 0; j < size1; ++j)
        (*this)(i,j) = val;

    return val;
  }

protected:
  T *v_;       // value
  int m_, n_;  // m X n dimensions of each matrix
  int Msize_;  // size of an individual matrix
  int size0_;   // number of matrices in the 1st dimension of the array
  int size1_;   // number of matrices in the 2nd dimension of the array
};

template< class T >
class MatrixD_Array2D : public MatrixDView_Array2D<T>
{
public:
  MatrixD_Array2D( const int m, const int n, const int size0, const int size1 ) :
    MatrixDView_Array2D<T>(new T[m*n*size0*size1], m, n, size0, size1)
  {
  }

  MatrixD_Array2D( const MatrixD_Array2D& A ) :
    MatrixDView_Array2D<T>(new T[A.m_*A.n_*A.size0_*A.size1_], A.m_, A.n_, A.size0_, A.size1_)
  {
    const int total = A.m_*A.n_*A.size0_*A.size1_;
    for (int i = 0; i < total; i++)
      v_[i] = A.v_[i];
  }

  T operator=( const T& val ) { return MatrixDView_Array2D<T>::operator=(val); }

  ~MatrixD_Array2D()
  {
    delete [] v_;
  }
protected:
  using MatrixDView_Array2D<T>::v_;
};

} //namespace SANS
} //namespace DLA

#endif // MATRIXD_ARRAY2D_H
