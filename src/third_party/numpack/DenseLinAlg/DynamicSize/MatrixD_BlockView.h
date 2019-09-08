// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_BLOCKVIEW_H
#define MATRIXD_BLOCKVIEW_H

#include <iostream>

#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"

#include "MatrixD_Type.h"
#include "numpack/DenseLinAlg/tools/Identity.h"
#include "tools/SANSTraitsPOD.h"


namespace numpack 
{
namespace DLA
{

//----------------------------------------------------------------------------//
// DenseBlockMatrixView:  Provieds a Block Matrix View of an array of Data
//
//----------------------------------------------------------------------------//

template< class T >
class DenseBlockMatrixView : public MatrixDType< DenseBlockMatrixView<T>, true >
{
private:
  // No default or copy constructor
  DenseBlockMatrixView() : v_(0), m_(0), n_(0), mB_(0), nB_(0), mS_(0), nS_(0), stride_(0), M_(v_, m_, n_, stride_) {}
  DenseBlockMatrixView( const DenseBlockMatrixView& z ) : v_(0), m_(0), n_(0), mB_(0), nB_(0), mS_(0), nS_(0), stride_(0), M_(v_, m_, n_, stride_) {}
  DenseBlockMatrixView& operator=( const DenseBlockMatrixView& ) = delete;

public:
  typedef T node_type;

  DenseBlockMatrixView( T v0[], const int m, const int n, const int mB, const int nB )
    : v_(v0), m_(m), n_(n), mB_(mB), nB_(nB), mS_(m / mB), nS_(n / nB), stride_(n), M_(v_, m_, n_, stride_)
  {
    SANS_ASSERT( m_ % mB_ == 0 );
    SANS_ASSERT( n_ % nB_ == 0 );
  }
  DenseBlockMatrixView( T v0[], const int m, const int n, const int stride, const int mB, const int nB )
    : v_(v0), m_(m), n_(n), mB_(mB), nB_(nB), mS_(m / mB), nS_(n / nB), stride_(stride), M_(v_, m_, n_, stride_)
  {
    SANS_ASSERT( m_ % mB_ == 0 );
    SANS_ASSERT( n_ % nB_ == 0 );
  }
  DenseBlockMatrixView( MatrixDView<T>& M, const int mB, const int nB )
    : v_(&M(0,0)), m_(M.m()), n_(M.n()), mB_(mB), nB_(nB), mS_(m_ / mB_), nS_(n_ / nB_), stride_(M.stride()), M_(v_, m_, n_, stride_)
  {
    SANS_ASSERT( m_ % mB_ == 0 );
    SANS_ASSERT( n_ % nB_ == 0 );
  }
  ~DenseBlockMatrixView() {}

  int size() const { return mB_*nB_; }
  int m() const { return mB_; }
  int n() const { return nB_; }
  int stride() const { return stride_; }

  // value accessor operators
  MatrixDView<T> operator()( const int i, const int j )
  {
    return MatrixDView<T>( v_ + i*mB_*stride_ + j*nB_, mS_, nS_, stride_ );
  }
  const MatrixDView<T> operator()( const int i, const int j ) const
  {
    return MatrixDView<T>( v_ + i*mB_*stride_ + j*nB_, mS_, nS_, stride_ );
  }

  // assignment
  DenseBlockMatrixView& operator=( const T& s ) { M_ = s; return *this; };
  DenseBlockMatrixView& operator=( const typename POD<T>::type& s ) { M_ = s; return *this; };
  DenseBlockMatrixView& operator=( const Identity& I ) { M_ = I; return *this; };

  // unary operators; no side effects
  const DenseBlockMatrixView& operator+() const { return *this; }

  // binary accumulation operators
  DenseBlockMatrixView& operator*=( const T& s ) { M_ *= s; return *this; };
  DenseBlockMatrixView& operator/=( const T& s ) { M_ /= s; return *this; };

  DenseBlockMatrixView& operator*=( const typename POD<T>::type& s ) { M_ *= s; return *this; };
  DenseBlockMatrixView& operator/=( const typename POD<T>::type& s ) { M_ /= s; return *this; };

  // input/output
  template<class> friend std::istream& operator>>( std::istream&, MatrixDView<T>& );
  template<class> friend std::ostream& operator<<( std::ostream&, const MatrixDView<T>& );

private:
  T *v_;         // value
  int m_, n_;    // m X n dimensions of the original matrix
  int mB_, nB_;  // mB X nB number of block matrices
  int mS_, nS_;  // mS X nS the size of the block matrices
  int stride_;   // Stride for when the memory is larger than m x n
  MatrixDView<T> M_;   // A matrix view just to the original memory to make some operations simpler
};

} //namespace DLA
} //namespace numpack 


#endif // MATRIXD_BLOCKVIEW_H
