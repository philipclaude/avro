// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_DIAG_H
#define MATRIXS_DIAG_H

#include <stdint.h> // uintptr_t

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
//#include <boost/mpl/assert.hpp>

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"
#include "MatrixS_Type.h"

namespace tinymat 
{
namespace DLA
{


//Represents a diagonal matrix
template<int M_, class T>
class MatrixSDiag : public MatrixSType< MatrixSDiag<M_,T>, true, true >
{

public:
  typedef T Ttype;

  static const int M = M_;
  static const int N = M_;

  explicit MatrixSDiag( const VectorS<M_,T>& V ) : V(V) {}

  //Operator to access the matrix values
  inline T operator()(const int i, const int j) const { return i == j ? V[i] : 0; }
  inline T diag(const int i) const { return V[i]; }

  // Lazy expression operations
  inline void value(const T& sgn, MatrixS<M,N,T>& res) const
  {
    res = 0;
    for (int i = 0; i < M; ++i)
      res(i,i) = sgn*V[i];
  }
  inline void plus(const T& sgn, MatrixS<M,N,T>& res) const
  {
    for (int i = 0; i < M; ++i)
      res(i,i) += sgn*V[i];
  }

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return V.ID(); }

private:
  const VectorS<M_,T>& V; //A vector representing the diagonal elements
};


//=============================================================================
//Operator to generate a datatype to represent a diagonal matrix
template< int M, class T >
inline MatrixSDiag<M,T>
diag(const VectorS<M,T>& V)
{
  return MatrixSDiag<M,T>( V );
}


//-----------------------------------------------------------------------------
// This specialization is for multiplication between a diagonal matrix a matrix, i.e.
// M = D*A
template<int M_, int N_, class T>
class OpMulS< MatrixSDiag<M_,T>, MatrixS<M_,N_,T> > : public MatrixSType< OpMulS< MatrixSDiag<M_,T>, MatrixS<M_,N_,T> >, true, true >
{
public:
  typedef T Ttype;

  static const int M = M_;
  static const int N = N_;

  OpMulS(const MatrixSDiag<M_,T>& dl, const MatrixS<M_,N_,T>& mr) :  dl(dl), mr(mr) {}

  inline void value(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) = sgn * dl.diag(i) * mr(i,j);
  }
  inline void plus(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) += sgn * dl.diag(i) * mr(i,j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  operator T() const
  {
    //BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
    //BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
    MatrixS<1,1,T> M(*this);
    return M(0,0);
  }

  //Means to access the left and right entries
  const MatrixSDiag<M_,T>& left()  const { return dl; }
  const MatrixS<M_,N_,T>&  right() const { return mr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSDiag<M_,T>& dl;
  const MatrixS<M_,N_,T>& mr;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between a matrix a diagonal matrix, i.e.
// M = A*D
template<int M_, int N_, class T>
class OpMulS< MatrixS<M_,N_,T>, MatrixSDiag<N_,T> > : public MatrixSType< OpMulS< MatrixS<M_,N_,T>, MatrixSDiag<N_,T> >, true, true >
{
public:
  typedef T Ttype;

  static const int M = M_;
  static const int N = N_;

  OpMulS(const MatrixS<M_,N_,T>& ml, const MatrixSDiag<N_,T>& dr ) :  ml(ml), dr(dr) {}

  inline void value(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) = sgn * ml(i,j) * dr.diag(j);
  }
  inline void plus(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) += sgn * ml(i,j) * dr.diag(j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  operator T() const
  {
    //BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
    //BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
    MatrixS<1,1,T> M(*this);
    return M(0,0);
  }

  //Means to access the left and right entries
  const MatrixS<M_,N_,T>&  left() const { return ml; }
  const MatrixSDiag<N_,T>& right()  const { return dr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixS<M_,N_,T>& ml;
  const MatrixSDiag<N_,T>& dr;
};


//-----------------------------------------------------------------------------
// This specialization is for multiplication between a diagonal matrix a transposed matrix, i.e.
// M = D*Transpose(A)
template<int M_, int N_, class T>
class OpMulS< MatrixSDiag<M_,T>, MatrixSTranspose<N_,M_,T> >
  : public MatrixSType< OpMulS< MatrixSDiag<M_,T>, MatrixSTranspose<N_,M_,T> >, true, true >
{
public:
  typedef T Ttype;

  static const int M = M_;
  static const int N = N_;

  OpMulS(const MatrixSDiag<M_,T>& dl, const MatrixSTranspose<N_,M_,T>& mr) :  dl(dl), mr(mr) {}

  inline void value(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) = sgn * dl.diag(i) * mr(i,j);
  }
  inline void plus(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) += sgn * dl.diag(i) * mr(i,j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  operator T() const
  {
    //BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
    //BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
    MatrixS<1,1,T> M(*this);
    return M(0,0);
  }

  //Means to access the left and right entries
  const MatrixSDiag<M_,T>& left()  const { return dl; }
  const MatrixSTranspose<N_,M_,T>&  right() const { return mr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSDiag<M_,T>& dl;
  const MatrixSTranspose<N_,M_,T>& mr;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between a matrix a diagonal matrix, i.e.
// M = Transpose(A)*D
template<int M_, int N_, class T>
class OpMulS< MatrixSTranspose<N_,M_,T>, MatrixSDiag<N_,T> >
  : public MatrixSType< OpMulS< MatrixSTranspose<N_,M_,T>, MatrixSDiag<N_,T> >, true, true >
{
public:
  typedef T Ttype;

  static const int M = M_;
  static const int N = N_;

  OpMulS(const MatrixSTranspose<N_,M_,T>& ml, const MatrixSDiag<N_,T>& dr ) :  ml(ml), dr(dr) {}

  inline void value(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) = sgn * ml(i,j) * dr.diag(j);
  }
  inline void plus(const T& sgn, MatrixS<M, N, T>& res) const
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        res(i,j) += sgn * ml(i,j) * dr.diag(j);
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  operator T() const
  {
    //BOOST_MPL_ASSERT_RELATION( M, ==, 1 );
    //BOOST_MPL_ASSERT_RELATION( N, ==, 1 );
    MatrixS<1,1,T> M(*this);
    return M(0,0);
  }

  //Means to access the left and right entries
  const MatrixSTranspose<N_,M_,T>&  left() const { return ml; }
  const MatrixSDiag<N_,T>& right()  const { return dr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSTranspose<N_,M_,T>& ml;
  const MatrixSDiag<N_,T>& dr;
};

} //namespace DLA
} //namespace tinymat 



#endif //MATRIXS_DIAG_H
