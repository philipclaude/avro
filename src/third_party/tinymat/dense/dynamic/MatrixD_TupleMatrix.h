// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXD_TUPLE_MATRIX_H
#define MATRIXD_TUPLE_MATRIX_H

#include <stdint.h> // uintptr_t

#include "tools/SANSException.h"
#include "tools/SANSTraitsPOD.h"
#include "tinymat/dense/tools/Identity.h"
#include "MatrixD_Type.h"

namespace tinymat 
{
namespace DLA
{

//Forward declaration
template< class ExprL, class ExprR >
class MatrixDExpressionTuple;


//This is a tuple designed specifically to be a tuple of MatrixD types. This allows
//the expression (M1, M2) which are m X n1 and m X n2 matrices be treated as an m X (n1+n2) matrix
template< class MatrixL >
class MatrixDTuple : public MatrixDType< MatrixDTuple<MatrixL>, true >
{

public:
  typedef typename MatrixL::node_type node_type;
  typedef node_type T;

  template<class, class> friend class MatrixDExpressionTuple;

  MatrixDTuple( MatrixL& L, MatrixDView<T>& R ) : L(L), R(R),
                                                  Lm_(L.m()), Ln_(L.n()),
                                                  Rm_(R.m()), Rn_(R.n()),
                                                  m_(Lm_), n_(Ln_ + Rn_)
  {
    //The tuple only works for matrices with the same number of rows.
    SANS_ASSERT( L.m() + R.m() );
  }

  //Operator to access the matrix values
  inline       node_type& operator()(const int i, const int j)       { return j < Ln_ ? L(i,j) : R(i, j-Ln_); }
  inline const node_type& operator()(const int i, const int j) const { return j < Ln_ ? L(i,j) : R(i, j-Ln_); }

  // Assignment operators
  MatrixDTuple& operator=( const MatrixDView<T>& M )
  {
    SANS_ASSERT( m_ == M.m() );
    SANS_ASSERT( n_ == M.n() );

    MatrixDView<T> LM( const_cast<T*>(&M(0,   0)), Lm_, Ln_, M.stride() );
    MatrixDView<T> RM( const_cast<T*>(&M(0, Ln_)), Rm_, Rn_, M.stride() );
    L = LM;
    R = RM;
    return *this;
  }

  MatrixDTuple& operator+=( const MatrixDView<T>& M )
  {
    SANS_ASSERT( m_ == M.m() );
    SANS_ASSERT( n_ == M.n() );

    MatrixDView<T> LM( const_cast<T*>(&M(0,   0)), Lm_, Ln_, M.stride() );
    MatrixDView<T> RM( const_cast<T*>(&M(0, Ln_)), Rm_, Rn_, M.stride() );
    L += LM;
    R += RM;
    return *this;
  }

  MatrixDTuple& operator-=( const MatrixDView<T>& M )
  {
    SANS_ASSERT( m_ == M.m() );
    SANS_ASSERT( n_ == M.n() );

    MatrixDView<T> LM( const_cast<T*>(&M(0,   0)), Lm_, Ln_, M.stride() );
    MatrixDView<T> RM( const_cast<T*>(&M(0, Ln_)), Rm_, Rn_, M.stride() );
    L -= LM;
    R -= RM;
    return *this;
  }

  //Scalar assignment operator
  MatrixDTuple& operator=( const T& s ) { L = s; R = s; return *this; }
  MatrixDTuple& operator=( const typename POD<T>::type& s )  { L = s; R = s; return *this; }

  //Identity operator to make the matrix typle Identity
  MatrixDTuple& operator=( const Identity& I )
  {
    L = 0;
    R = 0;
    for ( int i = 0; i < n_; i++ )
      (*this)(i,i) = I;
    return *this;
  }

  //Scalar compound assign
  MatrixDTuple& operator*=( const T& s ) { L *= s; R *= s; return *this; }
  MatrixDTuple& operator*=( const typename POD<T>::type& s )  { L *= s; R *= s; return *this; }

  //Lazy expression assignment operators
  template<class Expr, bool useRF> MatrixDTuple& operator= ( const MatrixDType<Expr, useRF>& );
  template<class Expr, bool useRF> MatrixDTuple& operator+=( const MatrixDType<Expr, useRF>& );
  template<class Expr, bool useRF> MatrixDTuple& operator-=( const MatrixDType<Expr, useRF>& );
  template<class Expr, bool useRF> MatrixDTuple& operator*=( const MatrixDType<Expr, useRF>& );

  // Row manipulations
  inline void swap_rows(const int i, const int j )
  {
    L.swap_rows(i,j);
    R.swap_rows(i,j);
  }
  inline void scale_row(const int i, const T& a, const int start = 0, int end = 0)
  {
    L.scale_row(i, a, start, end < Ln_ ? end : Ln_);
    if ( end > Ln_ || end == 0)
      R.scale_row(i, a, start, std::max(end-Ln_,0));
  }
  inline void axpy_rows(const int i, const int j, const T& a, const int start = 0, int end = 0)
  {
    L.axpy_rows(i, j, a, start, end < Ln_ ? end : Ln_);
    if ( end > Ln_ || end == 0)
      R.axpy_rows(i, j, a, start, std::max(end-Ln_,0));
  }

  // Lazy expression operations
  inline void value(const Real sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m_ == res.m() );
    SANS_ASSERT( n_ == res.n() );

    MatrixDView<T> resL(&res(0,  0), Lm_, Ln_, res.stride() );
    MatrixDView<T> resR(&res(0,Ln_), Rm_, Rn_, res.stride() );
    L.value(sgn, resL);
    R.value(sgn, resR);
  }
  inline void plus(const Real sgn, MatrixDView<T>& res) const
  {
    SANS_ASSERT( m_ == res.m() );
    SANS_ASSERT( n_ == res.n() );

    MatrixDView<T> resL(&res(0,  0), Lm_, Ln_, res.stride() );
    MatrixDView<T> resR(&res(0,Ln_), Rm_, Rn_, res.stride() );
    L.plus(sgn, resL);
    R.plus(sgn, resR);
  }

  //Lazy expressions where the LHS is a MatrixTyple, i.e. (a, b) = (c, d);
  // a and c must be the same dimensions, b and d must be the same dimensions.
  // c and d could be expressions
  template< class ResMatrixL >
  inline void value(const Real sgn, MatrixDTuple<ResMatrixL>& res) const
  {
    L.value(sgn, res.L);
    R.value(sgn, res.R);
  }
  template< class ResMatrixL >
  inline void plus(const Real sgn, MatrixDTuple<ResMatrixL>& res) const
  {
    L.plus(sgn, res.L);
    R.plus(sgn, res.R);
  }

  //Comma operator to recursively generate tuples, i.e. (a, b, c).
  inline MatrixDTuple< MatrixDTuple<MatrixL> >
  operator,( MatrixDView<T>& Matrix )
  {
    return MatrixDTuple< MatrixDTuple<MatrixL> >(*this, Matrix);
  }

  //This generates an expression tuple if a matrix expression follows
  //a matrix tuple, i.e. (a, b, c + d)
  template< class Expr, bool useRF >
  inline MatrixDExpressionTuple< MatrixDTuple<MatrixL>, Expr >
  operator,( const MatrixDType<Expr, useRF>& e )
  {
    return MatrixDExpressionTuple< MatrixDTuple<MatrixL>, Expr >(*this, e.cast());
  }

  int size() const { return m_*n_; }
  int m() const { return m_; }
  int n() const { return n_; }

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return L.ID(); }

  //Give MatrixD access so value/plus can access members
  friend class MatrixDView<T>;
private:
  MatrixL& L; //L could be a Matrix, or a MatrixTuple
  MatrixDView<T>& R; //By definition, R is always a matrix

  int Lm_, Ln_;   // m X n dimensions of the left matrix
  int Rm_, Rn_;   // m X n dimensions of the right matrix

  int m_, n_;   // m X n dimensions of the matrix
};


template< class MatrixL >
template< class Expr, bool useRF >
inline MatrixDTuple<MatrixL>&
MatrixDTuple<MatrixL>::operator=( const MatrixDType<Expr, useRF>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.value(1, *this);

  return *this;
}

template< class MatrixL >
template< class Expr, bool useRF >
inline MatrixDTuple<MatrixL>&
MatrixDTuple<MatrixL>::operator+=( const MatrixDType<Expr, useRF>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.plus(1, *this);

  return *this;
}

template< class MatrixL >
template< class Expr, bool useRF >
inline MatrixDTuple<MatrixL>&
MatrixDTuple<MatrixL>::operator-=( const MatrixDType<Expr, useRF>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );
  SANS_ASSERT( n_ == Tree.n() );

  Tree.plus(-1, *this);

  return *this;
}

template< class MatrixL >
template< class Expr, bool useRF >
inline MatrixDTuple<MatrixL>&
MatrixDTuple<MatrixL>::operator*=( const MatrixDType<Expr, useRF>& r )
{
  MatrixD<T> tmp(*this);
  *this = tmp*r;
  return *this;
}


//Overloaded comma operator to generate a matrix tuple with the expression (a, b)
template<class T>
inline MatrixDTuple< MatrixDView<T> >
operator,( const MatrixDView<T>& L, const MatrixDView<T>& R)
{
  typedef MatrixDView<T> MatrixType;
  return MatrixDTuple<MatrixType>( const_cast<MatrixType&>(L)
                                 , const_cast<MatrixType&>(R) );
}

} //namespace DLA
} //namespace tinymat 


#endif //MATRIXD_TUPLE_MATRIX_H
