// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMD_CLASS_H
#define MATRIXSYMD_CLASS_H

// Symmetric Matrix class with compile-time size

// NOTES:
// - indexing: zero-based
// - C-style matrix storage (row first)
// - Stores the Lower triangular matrix

#include <iostream>
#include <stdint.h> // uintptr_t

#include "MatrixD_Type.h"
#include "MatrixD_Add.h"
#include "MatrixD_Mul.h"
#include "MatrixD.h"
#include "MatrixD_Transpose.h"
#include "MatrixD_Diag.h"
//#include "MatrixSymD_InverseCholesky.h"
//#include "MatrixSymD_InverseLDLT.h"
#include "Eigen.h"

#include "numpack/DenseLinAlg/tools/Matrix_Util.h"
#include "numpack/DenseLinAlg/tools/Identity.h"

#include "tools/SANSTraitsPOD.h"
#include "tools/SANSException.h"

#include <vector>

namespace numpack
{
namespace DLA
{

//----------------------------------------------------------------------------//
// MatrixSymD<T>: symmetric matrix w/ MxM data elements, each of type T
//
// template parameters:
//   MxM          matrix dimension
//   T            matrix data type (e.g. double)
//
// Notes:
//  -implemented as lower triangular 1-D array with C-style (row-first) matrix storage
//----------------------------------------------------------------------------//

template <class T>
//class MatrixSymD : public MatrixDType< MatrixSymD<T>, true >
class MatrixSymD : public MatrixDType< MatrixSymD<T>, true >
{
public:

  typedef T node_type;
  typedef T Ttype;

  int m() const { return m_; }
  int n() const { return m_; }

  //Default constructor does not initialize the data. This makes it behave more like POD
  //and valgrind can catch any use of uninitialized data
  MatrixSymD() {}
  MatrixSymD( const int _m )
    : m_(_m)
  {
    allocate();
  }
  MatrixSymD( const MatrixSymD& m );
  MatrixSymD( const Identity& I ) { *this = I; };
  MatrixSymD( const T s , int m ); // missing 'explicit' allows MatrixSymD<N,Real> q = 0
  MatrixSymD( const T s[], int n );
  MatrixSymD( const std::initializer_list< std::initializer_list<T> >& s ) { *this = s; }
  MatrixSymD( const typename numpack::POD<T>::type& s , int m );
  MatrixSymD( const EigenSystemPair<T>& LE ) { *this = LE; }
  ~MatrixSymD() {}

#define INDEX(i,j) ( i > j  ?  i*(i+1)/2 + j  :  j*(j+1)/2 + i )
  // accessor operators
        T& operator()( int i, int j )       { return data[INDEX(i,j)]; }
  const T& operator()( int i, int j ) const { return data[INDEX(i,j)]; }
#undef INDEX

  void allocate()
  {
    M = m_;
    N = M;
    SIZE = m_*(m_+1)/2;
    data.resize(SIZE,T(0));
  }

  // assignment
  MatrixSymD& operator=( const MatrixSymD& m );
  MatrixSymD& operator=( const T& s );
  MatrixSymD& operator=( const typename numpack::POD<T>::type& s );
  MatrixSymD& operator=( const Identity& I );
  MatrixSymD& operator=( const std::initializer_list< std::initializer_list<T> >& s );
  MatrixSymD& operator=( const EigenSystemPair<T>& LE );

  // unary operators; no side effects
  const MatrixSymD& operator+() const;
  const MatrixSymD  operator-() const;

  // binary accumulation operators
  MatrixSymD& operator+=( const MatrixSymD& m );
  MatrixSymD& operator-=( const MatrixSymD& m );
  MatrixSymD& operator+=( const EigenSystemPair<T>& LE );
  MatrixSymD& operator-=( const EigenSystemPair<T>& LE );
  MatrixSymD& operator*=( const T& s );
  MatrixSymD& operator/=( const T& s );
  MatrixSymD& operator*=( const typename numpack::POD<T>::type& s );
  MatrixSymD& operator/=( const typename numpack::POD<T>::type& s );

  // Lazy expression with recursive functions assignment and binary accumulation, i.e.
  // MatrixSymD<T> C = A + B;
  // C += A + B;
  // C -= A + B;

  // Lazy expression element-wise assignment and binary accumulation
  template<class Expr>             MatrixSymD( const MatrixDType<Expr, true>& r ) { *this = r; }
  template<class Expr> MatrixSymD& operator= ( const MatrixDType<Expr, true>& );
  template<class Expr> MatrixSymD& operator+=( const MatrixDType<Expr, true>& );
  template<class Expr> MatrixSymD& operator-=( const MatrixDType<Expr, true>& );

  // Operators for A^T*A
  template<class T1>             MatrixSymD( const OpMulD< MatrixDTranspose<T1>, MatrixD<T1> >& tree ) { *this = tree; }
  template<class T1> MatrixSymD& operator= ( const OpMulD< MatrixDTranspose<T1>, MatrixD<T1> >& tree ) { return assign2(tree); }
  template<class T1> MatrixSymD& operator+=( const OpMulD< MatrixDTranspose<T1>, MatrixD<T1> >& tree ) { return addAssign2(tree); }
  template<class T1> MatrixSymD& operator-=( const OpMulD< MatrixDTranspose<T1>, MatrixD<T1> >& tree ) { return subAssign2(tree); }

  // Operators for A*A^T
  template<class T1>             MatrixSymD( const OpMulD< MatrixD<T1>, MatrixDTranspose<T1> >& tree ) { *this = tree; }
  template<class T1> MatrixSymD& operator= ( const OpMulD< MatrixD<T1>, MatrixDTranspose<T1> >& tree ) { return assign2(tree); }
  template<class T1> MatrixSymD& operator+=( const OpMulD< MatrixD<T1>, MatrixDTranspose<T1> >& tree ) { return addAssign2(tree); }
  template<class T1> MatrixSymD& operator-=( const OpMulD< MatrixD<T1>, MatrixDTranspose<T1> >& tree ) { return subAssign2(tree); }

  // Operators for S*S, where S is a symmetric matrix
  template<class T1>             MatrixSymD( const OpMulD< MatrixSymD<T1>, MatrixSymD<T1> >& tree ) { *this = tree; }
  template<class T1> MatrixSymD& operator= ( const OpMulD< MatrixSymD<T1>, MatrixSymD<T1> >& tree ) { return assign2(tree); }
  template<class T1> MatrixSymD& operator+=( const OpMulD< MatrixSymD<T1>, MatrixSymD<T1> >& tree ) { return addAssign2(tree); }
  template<class T1> MatrixSymD& operator-=( const OpMulD< MatrixSymD<T1>, MatrixSymD<T1> >& tree ) { return subAssign2(tree); }

  // Operators for A^T*D*A
  template<class T1, class T2>
              MatrixSymD( const OpMulD< OpMulD<MatrixDTranspose<T1>,MatrixDDiag<T2>>, MatrixD<T1> >& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymD& operator= ( const OpMulD< OpMulD<MatrixDTranspose<T1>,MatrixDDiag<T2>>, MatrixD<T1> >& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator+=( const OpMulD< OpMulD<MatrixDTranspose<T1>,MatrixDDiag<T2>>, MatrixD<T1> >& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator-=( const OpMulD< OpMulD<MatrixDTranspose<T1>,MatrixDDiag<T2>>, MatrixD<T1> >& tree ) { return subAssign3(tree); }

  // Operators for A*D*A^T
  template<class T1, class T2>
              MatrixSymD( const OpMulD< OpMulD<MatrixD<T1>,MatrixDDiag<T2>>, MatrixDTranspose<T1> >& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymD& operator= ( const OpMulD< OpMulD<MatrixD<T1>,MatrixDDiag<T2>>, MatrixDTranspose<T1> >& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator+=( const OpMulD< OpMulD<MatrixD<T1>,MatrixDDiag<T2>>, MatrixDTranspose<T1> >& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator-=( const OpMulD< OpMulD<MatrixD<T1>,MatrixDDiag<T2>>, MatrixDTranspose<T1> >& tree ) { return subAssign3(tree); }

  // Operators for A^T*S*A, where S is a symmetric matrix
  template<class T1, class T2>
              MatrixSymD( const OpMulD< OpMulD< MatrixDTranspose<T1>, MatrixSymD<T2>>, MatrixD<T1>>& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymD& operator= ( const OpMulD< OpMulD< MatrixDTranspose<T1>, MatrixSymD<T2>>, MatrixD<T1>>& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator+=( const OpMulD< OpMulD< MatrixDTranspose<T1>, MatrixSymD<T2>>, MatrixD<T1>>& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator-=( const OpMulD< OpMulD< MatrixDTranspose<T1>, MatrixSymD<T2>>, MatrixD<T1>>& tree ) { return subAssign3(tree); }

  // Operators for A*S*A^T, where S is a symmetric matrix
  template<class T1, class T2>
              MatrixSymD( const OpMulD< OpMulD< MatrixD<T1>, MatrixSymD<T2>>, MatrixDTranspose<T1>>& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymD& operator= ( const OpMulD< OpMulD< MatrixD<T1>, MatrixSymD<T2>>, MatrixDTranspose<T1>>& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator+=( const OpMulD< OpMulD< MatrixD<T1>, MatrixSymD<T2>>, MatrixDTranspose<T1>>& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator-=( const OpMulD< OpMulD< MatrixD<T1>, MatrixSymD<T2>>, MatrixDTranspose<T1>>& tree ) { return subAssign3(tree); }

  // Operators for S0*S1*S0, where S0, S1 are symmetric matrices
  template<class T1, class T2>
              MatrixSymD( const OpMulD< OpMulD< MatrixSymD<T1>, MatrixSymD<T2>>, MatrixSymD<T1>>& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymD& operator= ( const OpMulD< OpMulD< MatrixSymD<T1>, MatrixSymD<T2>>, MatrixSymD<T1>>& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator+=( const OpMulD< OpMulD< MatrixSymD<T1>, MatrixSymD<T2>>, MatrixSymD<T1>>& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymD& operator-=( const OpMulD< OpMulD< MatrixSymD<T1>, MatrixSymD<T2>>, MatrixSymD<T1>>& tree ) { return subAssign3(tree); }

  // Operators for Cholesky inverse
  /*
  template<class Expr>
              MatrixSymD( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse ) { *this = inverse; }
  template<class Expr>
  MatrixSymD& operator= ( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse );
  template<class Expr>
  MatrixSymD& operator+=( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse );
  template<class Expr>
  MatrixSymD& operator-=( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse );

  // Operators for Cholesky solve
  template<class InvExpr, class RExpr>
              MatrixSymD( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve ) { *this = solve; }
  template<class InvExpr, class RExpr>
  MatrixSymD& operator= ( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymD& operator+=( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymD& operator-=( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve );

  // Operators for LDL^T inverse
  template<class Expr>
              MatrixSymD( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse ) { *this = inverse; }
  template<class Expr>
  MatrixSymD& operator= ( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse );
  template<class Expr>
  MatrixSymD& operator+=( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse );
  template<class Expr>
  MatrixSymD& operator-=( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse );

  // Operators for LDL^T solve
  template<class InvExpr, class RExpr>
              MatrixSymD( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve ) { *this = solve; }
  template<class InvExpr, class RExpr>
  MatrixSymD& operator= ( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymD& operator+=( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymD& operator-=( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve );

  */

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

  //Row manipulations
  inline void swap_rows(const int i, const int j );
  template< class Ta >
  inline void scale_row(const int i, const Ta& a, const int start = 0, int end = 0);
  template< class Ta >
  inline void axpy_rows(const int i, const int j, const Ta& a, const int start = 0, int end = 0);
  inline int max_row_in_col(const int j, const int start = 0) const;

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return (uintptr_t)this; }

  //Lazy Expression Operators
  // Assign the value of this matrix to res
  template<class Tres>
  //pcaplan inline void value(const T& sgn, MatrixD<Tres>& res) const
  inline void value(const T& sgn, MatrixDView<Tres>& res) const
  {
    if ( ID() == res.ID() && sgn == 1 ) return; //Assigning my this to this, nothing to do.

    SANS_ASSERT( (uintptr_t)this != (uintptr_t)&res ); //The variable on the left of the assignment cannot also appear on the right

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < M; ++j)
        res(i,j) = sgn*(*this)(i,j);
  }

  // Add the value of this matrix to res
  template<class Tres>
  //pcaplan inline void plus(const T& sgn, MatrixD<Tres>& res) const
  inline void plus(const T& sgn, MatrixDView<Tres>& res) const
  {
    SANS_ASSERT( ID() != res.ID() ); //The variable on the left of the assignment cannot also appear on the right

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < M; ++j)
        res(i,j) += sgn*(*this)(i,j);
  }

  //Element-wise expression for symmetric matrices
  inline const T& value(const int& i) const { return data[i]; }
  inline       T& value(const int& i)       { return data[i]; }

  std::vector<T> data;

protected:

  //T data[SIZE];
  //static const int M = 2;//M_;
  //static const int N = 2;// M_;
  //static const int SIZE = M*(M+1)/2;

  int m_;
  int M;
  int N;
  int SIZE;

  // Helper functions for all the variations of assignment operators that involve the transpose
  template< class Expr >
  MatrixSymD& assign2( const Expr& Tree );
  template< class Expr >
  MatrixSymD& addAssign2( const Expr& Tree );
  template< class Expr >
  MatrixSymD& subAssign2( const Expr& Tree );

  template< class Expr >
  MatrixSymD& assign3( const Expr& Tree );
  template< class Expr >
  MatrixSymD& addAssign3( const Expr& Tree );
  template< class Expr >
  MatrixSymD& subAssign3( const Expr& Tree );
};


// constructors

template <class T>
inline
MatrixSymD<T>::MatrixSymD( const MatrixSymD<T>& m0 ) :
  m_(m0.m())
{
  allocate();
  for (int i = 0; i < SIZE; i++)
    data[i] = m0.data[i];
}

template <class T>
inline
MatrixSymD<T>::MatrixSymD( const T s[], int n ) :
  m_(n)
{
  allocate();
  SANS_ASSERT(n == SIZE);
  for (n = 0; n < SIZE; n++)
    data[n] = s[n];
}

template <class T>
inline
MatrixSymD<T>::MatrixSymD( const T s , int m  ) :
  m_(m)
{
  allocate();
  for (int n = 0; n < SIZE; n++)
    data[n] = s;
}

// needed for MatrixSymD< Surreal<Z>>(Real)
template <class T>
inline
MatrixSymD<T>::MatrixSymD( const typename numpack::POD<T>::type& s , int m ) :
  m_(m)
{
  allocate();
  for (int n = 0; n < SIZE; n++)
    data[n] = s;
}

// assignment

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const MatrixSymD& m0 )
{
  m_ = m0.m();
  allocate();
  if (this != &m0)
  {
    for (int i = 0; i < SIZE; i++)
      data[i] = m0.data[i];
  }
  return *this;
}

template<class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] = s;
  return *this;
}

// needed for  MatrixSymD< Surreal> q; q = 0;
template<class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] = s;
  return *this;
}

template<class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const Identity& I )
{
  *this = 0;
  for (int i = 0; i < M; i++)
    (*this)(i,i) = I;

  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const std::initializer_list< std::initializer_list<T> >& s )
{
  SANS_ASSERT(s.size() == (std::size_t)M);
  int n = 0;
  auto row = s.begin();
  for (std::size_t i = 0; i < M; ++i, row++)
  {
    SANS_ASSERT( i+1 == row->size() );
    auto col = row->begin();
    for (std::size_t j = 0; j < i+1; ++j, col++)
      data[n++] = *col;
  }

  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const EigenSystemPair<T>& LE )
{
  //Compute the symmetric matrix from the Eigen values and vectors
  *this = LE.E*diag(LE.L)*Transpose(LE.E);

  return *this;
}

// unary operators; no side effects

template <class T>
inline const MatrixSymD<T>&
MatrixSymD<T>::operator+() const
{
  return *this;
}

template <class T>
inline const MatrixSymD<T>
MatrixSymD<T>::operator-() const
{
  MatrixSymD<T> m(N);
  for (int i = 0; i < SIZE; i++)
    m.data[i] = -data[i];
  return m;
}


// binary accumulation operators

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const MatrixSymD& m )
{
  for (int i = 0; i < SIZE; i++)
    data[i] += m.data[i];
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const MatrixSymD& m )
{
  for (int i = 0; i < SIZE; i++)
    data[i] -= m.data[i];
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const EigenSystemPair<T>& LE )
{
  *this += MatrixSymD<T>(LE);
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const EigenSystemPair<T>& LE )
{
  *this -= MatrixSymD<T>(LE);
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator*=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] *= s;
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator/=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] /= s;
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator*=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] *= s;
  return *this;
}

template <class T>
inline MatrixSymD<T>&
MatrixSymD<T>::operator/=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] /= s;
  return *this;
}


//-----------------------------------------------------------------------------
// Lazy expression with recursive functions assignment and binary accumulation
// Theses functions are called from the specific interfaces involving Transpose.
template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::assign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) = tmp(i,j);

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::addAssign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) += tmp(i,j);

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::subAssign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) -= tmp(i,j);

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::assign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) = tmp(i,j);

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::addAssign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) += tmp(i,j);

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::subAssign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  //pcaplan SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  m_ = Tree.m();
  allocate();

  MatrixD<T> tmp(m_,m_);
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) -= tmp(i,j);

  return *this;
}


//-----------------------------------------------------------------------------
// Element-wise lazy expression assignment and binary accumulation
template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  #if 0
  for ( int i = 0; i < SIZE; i++ )
    data[i] = Tree.value(i);
  #else
  m_ = r.m();
  allocate();
  MatrixD<T> tmp(r.m(),r.m());
  Tree.value(1., tmp);
  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) = tmp(i,j);
  #endif

  return *this;
}

template <class T>
template< class Expr>
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  //BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  //BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  #if 0
  for ( int i = 0; i < SIZE; i++ )
    data[i] += Tree.value(i);
  #else
  m_ = r.m();
  allocate();
  MatrixD<T> tmp(r.m(),r.m());
  Tree.value(1., tmp);
  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) += tmp(i,j);
  #endif

  return *this;
}

template <class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const MatrixDType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  //BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  //BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  #if 0
  for ( int i = 0; i < SIZE; i++ )
    data[i] -= Tree.value(i);
  #else
  m_ = r.m();
  allocate();
  MatrixD<T> tmp(r.m(),r.m());
  Tree.value(1., tmp);
  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) -= tmp(i,j);
  #endif

  return *this;
}

/*

//-----------------------------------------------------------------------------
// Lazy expression of matrix inverse based on Cholesky decomposition
template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse )
{
  inverse.value(1., *this);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse )
{
  MatrixSymD<T> tmp;
  inverse.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const Fixed::MatrixInverse< MatrixSymDCholeskySolver, Expr >& inverse )
{
  MatrixSymD<T> tmp;
  inverse.value(1., tmp);
  *this -= tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve )
{
  solve.value(1., *this);

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve )
{
  MatrixSymD<T> tmp;
  solve.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const Fixed::MatrixSolve< MatrixSymDCholeskySolver, InvExpr, RExpr >& solve )
{
  MatrixSymD<T> tmp;
  solve.value(1., tmp);
  *this -= tmp;

  return *this;
}


//-----------------------------------------------------------------------------
// Lazy expression of matrix inverse based on LDL^T decomposition
template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse )
{
  inverse.value(1., *this);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse )
{
  MatrixSymD<T> tmp;
  inverse.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const Fixed::MatrixInverse< MatrixSymDLDLTSolver, Expr >& inverse )
{
  MatrixSymD<T> tmp;
  inverse.value(1., tmp);
  *this -= tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator=( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve )
{
  solve.value(1., *this);

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator+=( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve )
{
  MatrixSymD<T> tmp;
  solve.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymD<T>&
MatrixSymD<T>::operator-=( const Fixed::MatrixSolve< MatrixSymDLDLTSolver, InvExpr, RExpr >& solve )
{
  MatrixSymD<T> tmp;
  solve.value(1., tmp);
  *this -= tmp;

  return *this;
}

*/

//-----------------------------------------------------------------------------
// debug dump of private data
template <class T>
void
MatrixSymD<T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixSymD<" << M << "," << N << ",T>:" << std::endl;
#if 1
  out << indent << "  data = ";
  for (int n = 0; n < SIZE; n++)
    out << data[n] << " ";
  out << std::endl;
#else     // only works for class T with member function dump()
  for (int n = 0; n < SIZE; n++)
  {
    out << indent << "  data[" << n << "] = ";
    data[n].dump(indentSize);
  }
#endif
}

//-----------------------------------------------------------------------------
    //Row operations

//Swaps rows i and j
template<class T>
inline void
MatrixSymD<T>::swap_rows(const int i, const int j )
{
  if ( i == j ) return; //Nothing to do...
  MatrixUtil_Native<T,T>::swap(&(*this)(i,0), &(*this)(j,0), i+1);
}

//Scales the row i
template< class T >
template< class Ta >
inline void
MatrixSymD<T>::scale_row(const int i, const Ta& a, const int start, int end)
{
  if (end == 0) end = i+1;
  MatrixUtil_Native<T,Ta>::scal(&(*this)(i,0) + start, a, end - start);
}

//Performs y = a*x + y where x = row(i) and y = row(j)
template< class T >
template< class Ta >
inline void
MatrixSymD<T>::axpy_rows(const int i, const int j, const Ta& a, const int start, int end)
{
  SANS_DEVELOPER_EXCEPTION("MatrixSymD<T>::axpy_rows Not implemented yet");

  //if (end == 0) end = i+1;
  //if ( i > j )
  //  for (int m = 0; m < j+1; m++)
  //    (*this)(j,m) += a*(*this)(i,m);
  //else
  //  for (int m = 0; m < i+1; m++)
  //    (*this)(j,m) += a*(*this)(i,m);
}

//Finds the row with maximum value in a column
template<class T>
inline int
MatrixSymD<T>::max_row_in_col(const int j, const int start) const
{
  SANS_ASSERT(false); //Not implemented yet...
  return 0;
  //return MatrixUtil_Native<T>::max_row_in_col(&(*this)(0,j), M, N, start);
}

//- - - - - - - - - - - - - - non-member functions - - - - - - - - - - - - - -//

// I/O

template <class T>
std::ostream&
operator<<( std::ostream& out, const MatrixSymD<T>& m )
{
  const int M = m.m();
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
      out << m(i,j) << ", ";
    out << std::endl;
  }
  return out;
}

template <class T>
std::istream&
operator>>( std::istream &in, MatrixSymD<T>& m )
{
  const int M = m.m();
  for (int i = 0; i < M; i++)
    for (int j = 0; j < M; j++)
      in >> m(i,j);
  return in;
}


} //namespace DLA
} //namespace numpack

#endif // MATRIXSYMS_CLASS_H
