// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_2X2_H
#define MATRIXBLOCK_2X2_H

#include <initializer_list>

#include "BlockLinAlg_Mul.h"
#include "BlockLinAlg_Add.h"

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include "numpack/DenseLinAlg/tools/Identity.h"

//----------------------------------------------------------------------------//
// A class for generically storing 2x2 block matrices where the
// Matrices can be of mixed types
//
//----------------------------------------------------------------------------//

namespace numpack 
{

namespace BLA
{

template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
class MatrixBlock_2x2 : public BlockMatrixType< MatrixBlock_2x2<Matrix00, Matrix01,
                                                                Matrix10, Matrix11> >,
                        public BlockLinAlgType< MatrixBlock_2x2<Matrix00, Matrix01,
                                                                Matrix10, Matrix11> >
{
public:
  typedef Matrix00 M00;
  typedef Matrix01 M01;
  typedef Matrix10 M10;
  typedef Matrix11 M11;

  typedef MatrixBlock_2x2< typename Matrix00::size_type, typename Matrix01::size_type,
                           typename Matrix10::size_type, typename Matrix11::size_type > size_type;

  template<class U00, class U01,
           class U10, class U11>
  MatrixBlock_2x2( const std::initializer_list<U00>& s00, const std::initializer_list<U01>& s01,
                   const std::initializer_list<U10>& s10, const std::initializer_list<U11>& s11) :
    m00(s00), m01(s01),
    m10(s10), m11(s11) {}

  template<class U00, class U01,
           class U10, class U11>
  MatrixBlock_2x2( const std::initializer_list< std::initializer_list<U00> >& s00, const std::initializer_list< std::initializer_list<U01> >& s01,
                   const std::initializer_list< std::initializer_list<U10> >& s10, const std::initializer_list< std::initializer_list<U11> >& s11) :
    m00(s00), m01(s01),
    m10(s10), m11(s11) {}

  MatrixBlock_2x2( const typename Matrix00::size_type& s00, const typename Matrix01::size_type& s01,
                   const typename Matrix10::size_type& s10, const typename Matrix11::size_type& s11) :
    m00(s00), m01(s01),
    m10(s10), m11(s11) {}

  // Generic constructor.
  template<class A00, class A01,
           class A10, class A11>
  MatrixBlock_2x2( const A00& a00, const A01& a01,
                   const A10& a10, const A11& a11 ) :
    m00(a00), m01(a01),
    m10(a10), m11(a11) {}


  // Generic constructor. For example, non-zero patterns used to construct sparse matricies
  template<class A00, class A01,
           class A10, class A11>
  MatrixBlock_2x2( const MatrixBlock_2x2<A00, A01,
                                         A10, A11>& A ) :
    m00(A.m00), m01(A.m01),
    m10(A.m10), m11(A.m11) {}

  // Copy constructor
  MatrixBlock_2x2( const MatrixBlock_2x2& A ) :
    m00(A.m00), m01(A.m01),
    m10(A.m10), m11(A.m11) {}

  //Assignment operators
  MatrixBlock_2x2& operator=(const Real& v);

  MatrixBlock_2x2& operator=(const MatrixBlock_2x2& A);

  MatrixBlock_2x2& operator=( const DLA::Identity& I);

  template<class T>
  MatrixBlock_2x2& operator*=( const T& s ) { m00 *= s; m01 *= s; m10 *= s; m11 *= s; return *this; }

  // Lazy expression with recursive functions assignment and binary accumulation
  template<class Expr> MatrixBlock_2x2& operator= ( const BlockLinAlgType<Expr>& );
  template<class Expr> MatrixBlock_2x2& operator+=( const BlockLinAlgType<Expr>& );
  template<class Expr> MatrixBlock_2x2& operator-=( const BlockLinAlgType<Expr>& );
  template<class Expr> MatrixBlock_2x2& operator*=( const BlockLinAlgType<Expr>& );

  void value(const Real sgn, MatrixBlock_2x2& A) const;
  void plus(const Real sgn, MatrixBlock_2x2& A) const;

  template<class Ta>
  inline void scale_row0(const Ta& a);
  template<class Ta>
  inline void scale_row1(const Ta& a);

  template<class Ta>
  inline void axpy_rows01(const Ta& a);
  template<class Ta>
  inline void axpy_rows10(const Ta& a);

  // Lazy expression functions
  template<class Vector0, class Vector1>
  void mulVec_value(const VectorBlock_2<Vector0, Vector1>& x,
                    const Real sgn,
                    VectorBlock_2<Vector0, Vector1>& b) const;

  template<class Vector0, class Vector1>
  void mulVec_plus(const VectorBlock_2<Vector0, Vector1>& x,
                   const Real sgn,
                   VectorBlock_2<Vector0, Vector1>& b) const;

  template<class A00, class A01, class A10, class A11>
  void mulMat_value(const MatrixBlock_2x2<A00, A01, A10, A11>& M,
                    const Real sgn,
                    MatrixBlock_2x2<A00, A01, A10, A11>& res) const;

  template<class A00, class A01, class A10, class A11>
  void mulMat_plus(const MatrixBlock_2x2<A00, A01, A10, A11>& M,
                   const Real sgn,
                   MatrixBlock_2x2<A00, A01, A10, A11>& res) const;

  int m() const { return 2; }
  int n() const { return 2; }

  Matrix00 m00;
  Matrix01 m01;
  Matrix10 m10;
  Matrix11 m11;
};

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator=(const Real& v)
{
  m00 = v;
  m01 = v;
  m10 = v;
  m11 = v;

  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator=(const MatrixBlock_2x2& A)
{
  m00 = A.m00;
  m01 = A.m01;
  m10 = A.m10;
  m11 = A.m11;

  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator=(const DLA::Identity& I)
{
  m00 = I;
  m01 = 0;
  m10 = 0;
  m11 = I;

  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template< class Expr >
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator=(const BlockLinAlgType<Expr>& expr)
{
  const Expr& Tree = expr.cast();
  Tree.value(1, *this);
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template< class Expr >
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator+=(const BlockLinAlgType<Expr>& expr)
{
  const Expr& Tree = expr.cast();
  Tree.plus(1, *this);
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template< class Expr >
inline MatrixBlock_2x2<Matrix00, Matrix01, Matrix10, Matrix11>&
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::operator-=(const BlockLinAlgType<Expr>& expr)
{
  const Expr& Tree = expr.cast();
  Tree.plus(-1, *this);
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::value(const Real sgn, MatrixBlock_2x2& res) const
{
  res.m00 = sgn*m00;
  res.m01 = sgn*m01;
  res.m10 = sgn*m10;
  res.m11 = sgn*m11;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::plus(const Real sgn, MatrixBlock_2x2& res) const
{
  res.m00 += sgn*m00;
  res.m01 += sgn*m01;
  res.m10 += sgn*m10;
  res.m11 += sgn*m11;
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Vector0, class Vector1>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulVec_value(const VectorBlock_2<Vector0, Vector1>& x,
                                                    const Real sgn,
                                                    VectorBlock_2<Vector0, Vector1>& b) const
{
  b = 0;
  mulVec_plus(x, sgn, b);
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Vector0, class Vector1>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulVec_plus(const VectorBlock_2<Vector0, Vector1>& x,
                                                   const Real sgn,
                                                   VectorBlock_2<Vector0, Vector1>& b) const
{
  b.v0 += sgn*(m00*x.v0 + m01*x.v1);
  b.v1 += sgn*(m10*x.v0 + m11*x.v1);
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class A00, class A01, class A10, class A11>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulMat_value(const MatrixBlock_2x2<A00, A01, A10, A11>& M,
                                                    const Real sgn,
                                                    MatrixBlock_2x2<A00, A01, A10, A11>& res) const
{
  res = 0;
  mulMat_plus(M, sgn, res);
}

//=============================================================================
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class A00, class A01, class A10, class A11>
void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulMat_plus(const MatrixBlock_2x2<A00, A01, A10, A11>& M,
                                                   const Real sgn,
                                                   MatrixBlock_2x2<A00, A01, A10, A11>& res) const
{
  res.m00 += sgn*(m00*M.m00 + m01*M.m10);
  res.m01 += sgn*(m00*M.m01 + m01*M.m11);
  res.m10 += sgn*(m10*M.m00 + m11*M.m10);
  res.m11 += sgn*(m10*M.m01 + m11*M.m11);
}

//---------------------------------------------------------------------------//
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Ta>
inline void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::scale_row0(const Ta& a)
{
  Matrix00 tmp00 = a*m00;
  Matrix01 tmp01 = a*m01;
  m00 = tmp00;
  m01 = tmp01;
}
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Ta>
inline void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::scale_row1(const Ta& a)
{
  Matrix10 tmp10 = a*m10;
  Matrix11 tmp11 = a*m11;
  m10 = tmp10;
  m11 = tmp11;
}


//Performs y = a*x + y where x = row(0) and y = row(1)
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Ta>
inline void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::axpy_rows01(const Ta& a)
{
  m10 += a*m00;
  m11 += a*m01;
}

//Performs y = a*x + y where x = row(1) and y = row(0)
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Ta>
inline void
MatrixBlock_2x2< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::axpy_rows10(const Ta& a)
{
  m00 += a*m10;
  m01 += a*m11;
}

} //namespace BLA
} //namespace numpack 

#endif //MATRIXBLOCK_2X2_H
