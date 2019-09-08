// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXSYMS_CLASS_H
#define MATRIXSYMS_CLASS_H

// Symmetric Matrix class with compile-time size

// NOTES:
// - indexing: zero-based
// - C-style matrix storage (row first)
// - Stores the Lower triangular matrix

#include <iostream>
#include <stdint.h> // uintptr_t

#include "MatrixS_Type.h"
#include "MatrixS_Add.h"
#include "MatrixS_Mul.h"
#include "MatrixS.h"
#include "MatrixS_Transpose.h"
#include "MatrixS_Diag.h"
#include "MatrixSymS_InverseCholesky.h"
#include "MatrixSymS_InverseLDLT.h"
#include "Eigen.h"

#include "numpack/DenseLinAlg/tools/Matrix_Util.h"
#include "numpack/DenseLinAlg/tools/Identity.h"

#include "tools/SANSTraitsPOD.h"
#include "tools/SANSException.h"

namespace numpack 
{
namespace DLA
{

//----------------------------------------------------------------------------//
// MatrixSymS<M,T>: symmetric matrix w/ MxM data elements, each of type T
//
// template parameters:
//   MxM          matrix dimension
//   T            matrix data type (e.g. double)
//
// Notes:
//  -implemented as lower triangular 1-D array with C-style (row-first) matrix storage
//----------------------------------------------------------------------------//

template <int M_, class T>
class MatrixSymS : public MatrixSType< MatrixSymS<M_, T>, true, false >
{
public:
  typedef T Ttype;
  static const int M = M_;
  static const int N = M_;
  static const int SIZE = M*(M+1)/2;

  //Default constructor does not initialize the data. This makes it behave more like POD
  //and valgrind can catch any use of uninitialized data
  MatrixSymS() {}
  MatrixSymS( const MatrixSymS& m );
  MatrixSymS( const Identity& I ) { *this = I; };
  MatrixSymS( const T s ); // missing 'explicit' allows MatrixSymS<M,N,Real> q = 0
  MatrixSymS( const T s[], int n );
  MatrixSymS( const std::initializer_list< std::initializer_list<T> >& s ) { *this = s; }
  MatrixSymS( const typename numpack::POD<T>::type& s );
  MatrixSymS( const EigenSystemPair<M,T>& LE ) { *this = LE; }
  ~MatrixSymS() {}

#define INDEX(i,j) ( i > j  ?  i*(i+1)/2 + j  :  j*(j+1)/2 + i )
  // accessor operators
        T& operator()( int i, int j )       { return data[INDEX(i,j)]; }
  const T& operator()( int i, int j ) const { return data[INDEX(i,j)]; }
#undef INDEX

  // assignment
  MatrixSymS& operator=( const MatrixSymS& m );
  MatrixSymS& operator=( const T& s );
  MatrixSymS& operator=( const typename numpack::POD<T>::type& s );
  MatrixSymS& operator=( const Identity& I );
  MatrixSymS& operator=( const std::initializer_list< std::initializer_list<T> >& s );
  MatrixSymS& operator=( const EigenSystemPair<M,T>& LE );

  // unary operators; no side effects
  const MatrixSymS& operator+() const;
  const MatrixSymS  operator-() const;

  // binary accumulation operators
  MatrixSymS& operator+=( const MatrixSymS& m );
  MatrixSymS& operator-=( const MatrixSymS& m );
  MatrixSymS& operator+=( const EigenSystemPair<M,T>& LE );
  MatrixSymS& operator-=( const EigenSystemPair<M,T>& LE );
  MatrixSymS& operator*=( const T& s );
  MatrixSymS& operator/=( const T& s );
  MatrixSymS& operator*=( const typename numpack::POD<T>::type& s );
  MatrixSymS& operator/=( const typename numpack::POD<T>::type& s );

  // Lazy expression with recursive functions assignment and binary accumulation, i.e.
  // MatrixSymS<M,T> C = A + B;
  // C += A + B;
  // C -= A + B;

  // Lazy expression element-wise assignment and binary accumulation
  template<class Expr>             MatrixSymS( const MatrixSType<Expr, true, false>& r ) { *this = r; }
  template<class Expr> MatrixSymS& operator= ( const MatrixSType<Expr, true, false>& );
  template<class Expr> MatrixSymS& operator+=( const MatrixSType<Expr, true, false>& );
  template<class Expr> MatrixSymS& operator-=( const MatrixSType<Expr, true, false>& );

  // Operators for A^T*A
  template<int N_, class T1>             MatrixSymS( const OpMulS< MatrixSTranspose<N_,M,T1>, MatrixS<N_,M,T1> >& tree ) { *this = tree; }
  template<int N_, class T1> MatrixSymS& operator= ( const OpMulS< MatrixSTranspose<N_,M,T1>, MatrixS<N_,M,T1> >& tree ) { return assign2(tree); }
  template<int N_, class T1> MatrixSymS& operator+=( const OpMulS< MatrixSTranspose<N_,M,T1>, MatrixS<N_,M,T1> >& tree ) { return addAssign2(tree); }
  template<int N_, class T1> MatrixSymS& operator-=( const OpMulS< MatrixSTranspose<N_,M,T1>, MatrixS<N_,M,T1> >& tree ) { return subAssign2(tree); }

  // Operators for A*A^T
  template<int N_, class T1>             MatrixSymS( const OpMulS< MatrixS<M,N_,T1>, MatrixSTranspose<M,N_,T1> >& tree ) { *this = tree; }
  template<int N_, class T1> MatrixSymS& operator= ( const OpMulS< MatrixS<M,N_,T1>, MatrixSTranspose<M,N_,T1> >& tree ) { return assign2(tree); }
  template<int N_, class T1> MatrixSymS& operator+=( const OpMulS< MatrixS<M,N_,T1>, MatrixSTranspose<M,N_,T1> >& tree ) { return addAssign2(tree); }
  template<int N_, class T1> MatrixSymS& operator-=( const OpMulS< MatrixS<M,N_,T1>, MatrixSTranspose<M,N_,T1> >& tree ) { return subAssign2(tree); }

  // Operators for S*S, where S is a symmetric matrix
  template<class T1>             MatrixSymS( const OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T1> >& tree ) { *this = tree; }
  template<class T1> MatrixSymS& operator= ( const OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T1> >& tree ) { return assign2(tree); }
  template<class T1> MatrixSymS& operator+=( const OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T1> >& tree ) { return addAssign2(tree); }
  template<class T1> MatrixSymS& operator-=( const OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T1> >& tree ) { return subAssign2(tree); }

  // Operators for A^T*D*A
  template<int N_, class T1, class T2>
              MatrixSymS( const OpMulS< OpMulS<MatrixSTranspose<N_,M,T1>,MatrixSDiag<N_,T2>>, MatrixS<N_,M,T1> >& tree ) { *this = tree; }
  template<int N_, class T1, class T2>
  MatrixSymS& operator= ( const OpMulS< OpMulS<MatrixSTranspose<N_,M,T1>,MatrixSDiag<N_,T2>>, MatrixS<N_,M,T1> >& tree ) { return assign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator+=( const OpMulS< OpMulS<MatrixSTranspose<N_,M,T1>,MatrixSDiag<N_,T2>>, MatrixS<N_,M,T1> >& tree ) { return addAssign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator-=( const OpMulS< OpMulS<MatrixSTranspose<N_,M,T1>,MatrixSDiag<N_,T2>>, MatrixS<N_,M,T1> >& tree ) { return subAssign3(tree); }

  // Operators for A*D*A^T
  template<int N_, class T1, class T2>
              MatrixSymS( const OpMulS< OpMulS<MatrixS<M,N_,T1>,MatrixSDiag<N_,T2>>, MatrixSTranspose<M,N_,T1> >& tree ) { *this = tree; }
  template<int N_, class T1, class T2>
  MatrixSymS& operator= ( const OpMulS< OpMulS<MatrixS<M,N_,T1>,MatrixSDiag<N_,T2>>, MatrixSTranspose<M,N_,T1> >& tree ) { return assign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator+=( const OpMulS< OpMulS<MatrixS<M,N_,T1>,MatrixSDiag<N_,T2>>, MatrixSTranspose<M,N_,T1> >& tree ) { return addAssign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator-=( const OpMulS< OpMulS<MatrixS<M,N_,T1>,MatrixSDiag<N_,T2>>, MatrixSTranspose<M,N_,T1> >& tree ) { return subAssign3(tree); }

  // Operators for A^T*S*A, where S is a symmetric matrix
  template<int N_, class T1, class T2>
              MatrixSymS( const OpMulS< OpMulS< MatrixSTranspose<N_,M,T1>, MatrixSymS<N_,T2>>, MatrixS<N_,M,T1>>& tree ) { *this = tree; }
  template<int N_, class T1, class T2>
  MatrixSymS& operator= ( const OpMulS< OpMulS< MatrixSTranspose<N_,M,T1>, MatrixSymS<N_,T2>>, MatrixS<N_,M,T1>>& tree ) { return assign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator+=( const OpMulS< OpMulS< MatrixSTranspose<N_,M,T1>, MatrixSymS<N_,T2>>, MatrixS<N_,M,T1>>& tree ) { return addAssign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator-=( const OpMulS< OpMulS< MatrixSTranspose<N_,M,T1>, MatrixSymS<N_,T2>>, MatrixS<N_,M,T1>>& tree ) { return subAssign3(tree); }

  // Operators for A*S*A^T, where S is a symmetric matrix
  template<int N_, class T1, class T2>
              MatrixSymS( const OpMulS< OpMulS< MatrixS<M,N_,T1>, MatrixSymS<N_,T2>>, MatrixSTranspose<M,N_,T1>>& tree ) { *this = tree; }
  template<int N_, class T1, class T2>
  MatrixSymS& operator= ( const OpMulS< OpMulS< MatrixS<M,N_,T1>, MatrixSymS<N_,T2>>, MatrixSTranspose<M,N_,T1>>& tree ) { return assign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator+=( const OpMulS< OpMulS< MatrixS<M,N_,T1>, MatrixSymS<N_,T2>>, MatrixSTranspose<M,N_,T1>>& tree ) { return addAssign3(tree); }
  template<int N_, class T1, class T2>
  MatrixSymS& operator-=( const OpMulS< OpMulS< MatrixS<M,N_,T1>, MatrixSymS<N_,T2>>, MatrixSTranspose<M,N_,T1>>& tree ) { return subAssign3(tree); }

  // Operators for S0*S1*S0, where S0, S1 are symmetric matrices
  template<class T1, class T2>
              MatrixSymS( const OpMulS< OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T2>>, MatrixSymS<M,T1>>& tree ) { *this = tree; }
  template<class T1, class T2>
  MatrixSymS& operator= ( const OpMulS< OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T2>>, MatrixSymS<M,T1>>& tree ) { return assign3(tree); }
  template<class T1, class T2>
  MatrixSymS& operator+=( const OpMulS< OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T2>>, MatrixSymS<M,T1>>& tree ) { return addAssign3(tree); }
  template<class T1, class T2>
  MatrixSymS& operator-=( const OpMulS< OpMulS< MatrixSymS<M,T1>, MatrixSymS<M,T2>>, MatrixSymS<M,T1>>& tree ) { return subAssign3(tree); }

  // Operators for Cholesky inverse
  template<class Expr>
              MatrixSymS( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse ) { *this = inverse; }
  template<class Expr>
  MatrixSymS& operator= ( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse );
  template<class Expr>
  MatrixSymS& operator+=( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse );
  template<class Expr>
  MatrixSymS& operator-=( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse );

  // Operators for Cholesky solve
  template<class InvExpr, class RExpr>
              MatrixSymS( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve ) { *this = solve; }
  template<class InvExpr, class RExpr>
  MatrixSymS& operator= ( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymS& operator+=( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymS& operator-=( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve );

  // Operators for LDL^T inverse
  template<class Expr>
              MatrixSymS( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse ) { *this = inverse; }
  template<class Expr>
  MatrixSymS& operator= ( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse );
  template<class Expr>
  MatrixSymS& operator+=( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse );
  template<class Expr>
  MatrixSymS& operator-=( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse );

  // Operators for LDL^T solve
  template<class InvExpr, class RExpr>
              MatrixSymS( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve ) { *this = solve; }
  template<class InvExpr, class RExpr>
  MatrixSymS& operator= ( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymS& operator+=( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve );
  template<class InvExpr, class RExpr>
  MatrixSymS& operator-=( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve );

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
  inline void value(const T& sgn, MatrixS<M,M,Tres>& res) const
  {
    if ( ID() == res.ID() && sgn == 1 ) return; //Assigning my this to this, nothing to do.

    SANS_ASSERT( (uintptr_t)this != (uintptr_t)&res ); //The variable on the left of the assignment cannot also appear on the right

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < M; ++j)
        res(i,j) = sgn*(*this)(i,j);
  }

  // Add the value of this matrix to res
  template<class Tres>
  inline void plus(const T& sgn, MatrixS<M,M,Tres>& res) const
  {
    SANS_ASSERT( ID() != res.ID() ); //The variable on the left of the assignment cannot also appear on the right

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < M; ++j)
        res(i,j) += sgn*(*this)(i,j);
  }

  //Element-wise expression for symmetric matrices
  inline const T& value(const int& i) const { return data[i]; }
  inline       T& value(const int& i)       { return data[i]; }

protected:
  T data[SIZE];

  // Helper functions for all the variations of assignment operators that involve the transpose
  template< class Expr >
  MatrixSymS& assign2( const Expr& Tree );
  template< class Expr >
  MatrixSymS& addAssign2( const Expr& Tree );
  template< class Expr >
  MatrixSymS& subAssign2( const Expr& Tree );

  template< class Expr >
  MatrixSymS& assign3( const Expr& Tree );
  template< class Expr >
  MatrixSymS& addAssign3( const Expr& Tree );
  template< class Expr >
  MatrixSymS& subAssign3( const Expr& Tree );
};


// constructors

template <int M, class T>
inline
MatrixSymS<M,T>::MatrixSymS( const MatrixSymS& m )
{
  for (int i = 0; i < SIZE; i++)
    data[i] = m.data[i];
}

template <int M, class T>
inline
MatrixSymS<M,T>::MatrixSymS( const T s[], int n )
{
  SANS_ASSERT(n == SIZE);
  for (n = 0; n < SIZE; n++)
    data[n] = s[n];
}

template <int M, class T>
inline
MatrixSymS<M,T>::MatrixSymS( const T s )
{
  for (int n = 0; n < SIZE; n++)
    data[n] = s;
}

// needed for MatrixSymS<M, Surreal<Z>>(Real)
template <int M, class T>
inline
MatrixSymS<M,T>::MatrixSymS( const typename numpack::POD<T>::type& s )
{
  for (int n = 0; n < SIZE; n++)
    data[n] = s;
}

// assignment

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const MatrixSymS& m )
{
  if (this != &m)
  {
    for (int i = 0; i < SIZE; i++)
      data[i] = m.data[i];
  }
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] = s;
  return *this;
}

// needed for  MatrixSymS<M, Surreal> q; q = 0;
template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] = s;
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const Identity& I )
{
  *this = 0;
  for (int i = 0; i < M; i++)
    (*this)(i,i) = I;

  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const std::initializer_list< std::initializer_list<T> >& s )
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

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const EigenSystemPair<M,T>& LE )
{
  //Compute the symmetric matrix from the Eigen values and vectors
  *this = LE.E*diag(LE.L)*Transpose(LE.E);

  return *this;
}

// unary operators; no side effects

template <int M, class T>
inline const MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+() const
{
  return *this;
}

template <int M, class T>
inline const MatrixSymS<M,T>
MatrixSymS<M,T>::operator-() const
{
  MatrixSymS<M,T> m;
  for (int i = 0; i < SIZE; i++)
    m.data[i] = -data[i];
  return m;
}


// binary accumulation operators

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const MatrixSymS& m )
{
  for (int i = 0; i < SIZE; i++)
    data[i] += m.data[i];
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const MatrixSymS& m )
{
  for (int i = 0; i < SIZE; i++)
    data[i] -= m.data[i];
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const EigenSystemPair<M,T>& LE )
{
  *this += MatrixSymS<M,T>(LE);
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const EigenSystemPair<M,T>& LE )
{
  *this -= MatrixSymS<M,T>(LE);
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator*=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] *= s;
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator/=( const T& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] /= s;
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator*=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] *= s;
  return *this;
}

template <int M, class T>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator/=( const typename numpack::POD<T>::type& s )
{
  for (int i = 0; i < SIZE; i++)
    data[i] /= s;
  return *this;
}


//-----------------------------------------------------------------------------
// Lazy expression with recursive functions assignment and binary accumulation
// Theses functions are called from the specific interfaces involving Transpose.
template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::assign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) = tmp(i,j);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::addAssign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) += tmp(i,j);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::subAssign2( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) -= tmp(i,j);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::assign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) = tmp(i,j);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::addAssign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) += tmp(i,j);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::subAssign3( const Expr& Tree )
{
  //Make sure the expression is of the form A^T*M*A which creates a symmetric matrix
  SANS_ASSERT( Tree.left().left().ID() == Tree.right().ID() );

  MatrixS<M,M,T> tmp;
  Tree.value(1., tmp);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < i+1; j++)
      (*this)(i,j) -= tmp(i,j);

  return *this;
}


//-----------------------------------------------------------------------------
// Element-wise lazy expression assignment and binary accumulation

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const MatrixSType<Expr, true, false>& r )
{
  const Expr& Tree = r.cast();

  //BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  //BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for ( int i = 0; i < SIZE; i++ )
    data[i] = Tree.value(i);

  return *this;
}

template <int M, class T>
template< class Expr>
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const MatrixSType<Expr, true, false>& r )
{
  const Expr& Tree = r.cast();

  //BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  //BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for ( int i = 0; i < SIZE; i++ )
    data[i] += Tree.value(i);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const MatrixSType<Expr, true, false>& r )
{
  const Expr& Tree = r.cast();

  //BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  //BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for ( int i = 0; i < SIZE; i++ )
    data[i] -= Tree.value(i);

  return *this;
}

//-----------------------------------------------------------------------------
// Lazy expression of matrix inverse based on Cholesky decomposition
template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse )
{
  inverse.value(1., *this);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse )
{
  MatrixSymS<M,T> tmp;
  inverse.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const Fixed::MatrixInverse< MatrixSymSCholeskySolver, Expr >& inverse )
{
  MatrixSymS<M,T> tmp;
  inverse.value(1., tmp);
  *this -= tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve )
{
  solve.value(1., *this);

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve )
{
  MatrixSymS<M,T> tmp;
  solve.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const Fixed::MatrixSolve< MatrixSymSCholeskySolver, InvExpr, RExpr >& solve )
{
  MatrixSymS<M,T> tmp;
  solve.value(1., tmp);
  *this -= tmp;

  return *this;
}


//-----------------------------------------------------------------------------
// Lazy expression of matrix inverse based on LDL^T decomposition
template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse )
{
  inverse.value(1., *this);

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse )
{
  MatrixSymS<M,T> tmp;
  inverse.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class Expr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const Fixed::MatrixInverse< MatrixSymSLDLTSolver, Expr >& inverse )
{
  MatrixSymS<M,T> tmp;
  inverse.value(1., tmp);
  *this -= tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator=( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve )
{
  solve.value(1., *this);

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator+=( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve )
{
  MatrixSymS<M,T> tmp;
  solve.value(1., tmp);
  *this += tmp;

  return *this;
}

template <int M, class T>
template< class InvExpr, class RExpr >
inline MatrixSymS<M,T>&
MatrixSymS<M,T>::operator-=( const Fixed::MatrixSolve< MatrixSymSLDLTSolver, InvExpr, RExpr >& solve )
{
  MatrixSymS<M,T> tmp;
  solve.value(1., tmp);
  *this -= tmp;

  return *this;
}

//-----------------------------------------------------------------------------
// debug dump of private data
template <int M, class T>
void
MatrixSymS<M,T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixS<" << M << "," << N << ",T>:" << std::endl;
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
template< int M, class T >
inline void
MatrixSymS<M,T>::swap_rows(const int i, const int j )
{
  if ( i == j ) return; //Nothing to do...
  MatrixUtil_Native<T,T>::swap(&(*this)(i,0), &(*this)(j,0), i+1);
}

//Scales the row i
template< int M, class T >
template< class Ta >
inline void
MatrixSymS<M,T>::scale_row(const int i, const Ta& a, const int start, int end)
{
  if (end == 0) end = i+1;
  MatrixUtil_Native<T,Ta>::scal(&(*this)(i,0) + start, a, end - start);
}

//Performs y = a*x + y where x = row(i) and y = row(j)
template< int M, class T >
template< class Ta >
inline void
MatrixSymS<M,T>::axpy_rows(const int i, const int j, const Ta& a, const int start, int end)
{
  SANS_DEVELOPER_EXCEPTION("MatrixSymS<M,T>::axpy_rows Not implemented yet");

  //if (end == 0) end = i+1;
  //if ( i > j )
  //  for (int m = 0; m < j+1; m++)
  //    (*this)(j,m) += a*(*this)(i,m);
  //else
  //  for (int m = 0; m < i+1; m++)
  //    (*this)(j,m) += a*(*this)(i,m);
}

//Finds the row with maximum value in a column
template< int M, class T >
inline int
MatrixSymS<M,T>::max_row_in_col(const int j, const int start) const
{
  SANS_ASSERT(false); //Not implemented yet...
  return 0;
  //return MatrixUtil_Native<T>::max_row_in_col(&(*this)(0,j), M, N, start);
}

//- - - - - - - - - - - - - - non-member functions - - - - - - - - - - - - - -//

// I/O

template <int M, class T>
std::ostream&
operator<<( std::ostream& out, const MatrixSymS<M,T>& m )
{
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
      out << m(i,j) << ", ";
    out << std::endl;
  }
  return out;
}

template <int M, class T>
std::istream&
operator>>( std::istream &in, MatrixSymS<M,T>& m )
{
  for (int i = 0; i < M; i++)
    for (int j = 0; j < M; j++)
      in >> m(i,j);
  return in;
}


} //namespace DLA
} //namespace numpack 

#endif // MATRIXSYMS_CLASS_H
