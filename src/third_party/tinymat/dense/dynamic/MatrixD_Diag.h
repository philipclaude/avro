// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_DIAG_H
#define MATRIXD_DIAG_H

#include <stdint.h> // uintptr_t

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixD_Type.h"

namespace tinymat 
{
namespace DLA
{


//Represents a diagonal matrix
template<class T>
class MatrixDDiag : public MatrixDType< MatrixDDiag<T>, true >
{

public:
  typedef T node_type;

  explicit MatrixDDiag( const VectorDView<T>& V ) : V_(V) {}

  //Operator to access the matrix values
  inline T operator()(const int i, const int j) const { return i == j ? V_[i] : 0; }
  inline T diag(const int i) const { return V_[i]; }

  // Lazy expression operations
  inline void value(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    res = 0;
    for (int i = 0; i < m; ++i)
      res(i,i) = sgn*V_[i];
  }
  inline void plus(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    for (int i = 0; i < m; ++i)
      res(i,i) += sgn*V_[i];
  }

  int size() const { return m()*n(); }
  int m() const { return V_.m(); }
  int n() const { return V_.m(); } //Diagonal square matrix, so n = V_.m
  int stride() const { return V_.stride(); } //Memory stride does not change

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return V_.ID(); }

private:
  const VectorDView<T>& V_; //A vector representing the diagonal elements
};


//=============================================================================
//Operator to generate a datatype to represent a diagonal matrix
template< class T >
inline MatrixDDiag<T>
diag(const VectorDView<T>& V)
{
  return MatrixDDiag<T>( V );
}


//-----------------------------------------------------------------------------
// This specialization is for multiplication between a diagonal matrix a matrix, i.e.
// M = D*A
template<class T>
class OpMulD< MatrixDDiag<T>, MatrixDView<T> > : public MatrixDType< OpMulD< MatrixDDiag<T>, MatrixDView<T> >, true >
{
public:
  typedef T node_type;

  OpMulD(const MatrixDDiag<T>& dl, const MatrixDView<T>& mr) :  dl_(dl), mr_(mr) {}

  inline void value(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    const int n = res.n();

    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        res(i,j) = sgn * dl_.diag(i) * mr_(i,j);
  }
  inline void plus(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    const int n = res.n();

    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        res(i,j) += sgn * dl_.diag(i) * mr_(i,j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<T> MR( e );
    MatrixD<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<T> MR( e );
    MatrixD<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

//  operator T() const
//  {
//    BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
//    BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
//    MatrixD<1,1,T> M(*this);
//    return M(0,0);
//  }

  int size() const { return m()*n(); }
  int m() const { return dl_.m(); }
  int n() const { return mr_.n(); }
  int stride() const { return mr_.stride(); } //Memory stride does not change

  //Means to access the left and right entries
  //const MatrixDDiag<M_,T>& left()  const { return dl; }
  //const MatrixD<M_,N_,T>&  right() const { return mr; }

  inline const OpMulD&
  operator+() const { return *this; }
private:
  const MatrixDDiag<T>& dl_;
  const MatrixDView<T>& mr_;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between a diagonal matrix a matrix, i.e.
// M = A*D
template<class T>
class OpMulD< MatrixDView<T>, MatrixDDiag<T> > : public MatrixDType< OpMulD< MatrixDView<T>, MatrixDDiag<T> >, true >
{
public:
  typedef T node_type;

  OpMulD(const MatrixDView<T>& ml, const MatrixDDiag<T>& dr) : ml_(ml), dr_(dr) {}

  inline void value(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    const int n = res.n();

    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        res(i,j) = sgn * ml_(i,j) * dr_.diag(j);
  }
  inline void plus(const T& sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m() == res.m() );
    SANS_ASSERT( n() == res.n() );

    const int m = res.m();
    const int n = res.n();

    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        res(i,j) += sgn * ml_(i,j) * dr_.diag(j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<T> MR( e );
    MatrixD<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixDTuple<MatrixL>& tuple) const
  {
    MatrixD<T> MR( e );
    MatrixD<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

//  operator T() const
//  {
//    BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
//    BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
//    MatrixD<1,1,T> M(*this);
//    return M(0,0);
//  }

  int size() const { return m()*n(); }
  int m() const { return ml_.m(); }
  int n() const { return dr_.n(); }
  int stride() const { return ml_.stride(); } //Memory stride does not change

  //Means to access the left and right entries
  //const MatrixDDiag<M_,T>& left()  const { return dl; }
  //const MatrixD<M_,N_,T>&  right() const { return mr; }

  inline const OpMulD&
  operator+() const { return *this; }
private:
  const MatrixDView<T>& ml_;
  const MatrixDDiag<T>& dr_;
};

} //namespace DLA
} //namespace tinymat 



#endif //MATRIXS_DIAG_H
