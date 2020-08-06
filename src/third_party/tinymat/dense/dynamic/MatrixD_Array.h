// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_ARRAY_H
#define MATRIXD_ARRAY_H

#include "tools/SANSnumerics.h"
#include "MatrixD.h"

namespace tinymat 
{
namespace DLA
{

//-----------------------------------------------------------------------------
//
// MatrixD_Array provides an array of equally sized matrices.
// The implementation uses continuous memory
//
//-----------------------------------------------------------------------------

template< class T >
class MatrixDView_Array
{
public:
  MatrixDView_Array( T v0[], const int m, const int n, const int size ) : v_(v0), m_(m), n_(n), Msize_(m*n), size_(size)
  {
    SANS_ASSERT( m_ >= 0 );
    SANS_ASSERT( n_ >= 0 );
  }
  MatrixDView_Array( const MatrixDView_Array& A ) : v_(A.v_), m_(A.m_), n_(A.n_), Msize_(A.Msize_), size_(A.size_) {}
  MatrixDView_Array& operator=( const MatrixDView_Array& ) = delete;

  int size() const { return size_; }
  int m() const { return m_; }
  int n() const { return n_; }

  // value accessor operators
        MatrixDView<T> operator[]( const int i )       { return MatrixDView<T>(v_ + i*Msize_, m_, n_); }
  const MatrixDView<T> operator[]( const int i ) const { return MatrixDView<T>(v_ + i*Msize_, m_, n_); }

  // element wise accessors
        T& value(const int i)       { return v_[i]; }
  const T& value(const int i) const { return v_[i]; }

  T operator=( const T& val )
  {
    for (int i = 0; i < size_; ++i)
      (*this)[i] = val;

    return val;
  }

protected:
  T *v_;       // value
  int m_, n_;  // m X n dimensions of each matrix
  int Msize_;  // size of an individual matrix
  int size_;   // number of matrices in the array
};

template< class T >
class MatrixD_Array : public MatrixDView_Array<T>
{
public:
  MatrixD_Array( const int m, const int n, const int size ) : MatrixDView_Array<T>(new T[m*n*size], m, n, size)
  {
  }

  MatrixD_Array( const MatrixD_Array& A ) : MatrixDView_Array<T>(new T[A.m_*A.n_*A.size_], A.m_, A.n_, A.size_)
  {
    const int total = A.m_*A.n_*A.size_;
    for (int i = 0; i < total; i++)
      v_[i] = A.v_[i];
  }

  T operator=( const T& val ) { return MatrixDView_Array<T>::operator=(val); }

  ~MatrixD_Array()
  {
    delete [] v_;
  }
protected:
  using MatrixDView_Array<T>::v_;
};

} //namespace tinymat 
} //namespace DLA

#endif // MATRIXD_ARRAY_H
