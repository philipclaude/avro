// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXBLOCK_3X3_H
#define MATRIXBLOCK_3X3_H

#include <initializer_list>

#include "block_Mul.h"
#include "block_Add.h"

#include "tools/SANSnumerics.h"

#include "numpack/dense/tools/Identity.h"

//----------------------------------------------------------------------------//
// A class for generically storing 3x3 block matrices where the
// Matrices can be of mixed types
//
//----------------------------------------------------------------------------//

namespace numpack 
{

namespace BLA
{

template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
class MatrixBlock_3x3 : public BlockMatrixType< MatrixBlock_3x3<Matrix00, Matrix01, Matrix02,
                                                                Matrix10, Matrix11, Matrix12,
                                                                Matrix20, Matrix21, Matrix22> >,
                        public blockType< MatrixBlock_3x3<Matrix00, Matrix01, Matrix02,
                                                                Matrix10, Matrix11, Matrix12,
                                                                Matrix20, Matrix21, Matrix22> >
{
public:

  typedef MatrixBlock_3x3< typename Matrix00::size_type, typename Matrix01::size_type, typename Matrix02::size_type,
                           typename Matrix10::size_type, typename Matrix11::size_type, typename Matrix12::size_type,
                           typename Matrix20::size_type, typename Matrix21::size_type, typename Matrix22::size_type> size_type;


  template<class U00, class U01, class U02,
           class U10, class U11, class U12,
           class U20, class U21, class U22>
  MatrixBlock_3x3( const std::initializer_list<U00>& s00, const std::initializer_list<U01>& s01, const std::initializer_list<U02>& s02,
                   const std::initializer_list<U10>& s10, const std::initializer_list<U11>& s11, const std::initializer_list<U12>& s12,
                   const std::initializer_list<U20>& s20, const std::initializer_list<U21>& s21, const std::initializer_list<U22>& s22) :
    m00(s00), m01(s01), m02(s02),
    m10(s10), m11(s11), m12(s12),
    m20(s20), m21(s21), m22(s22) {}

  template<class U00, class U01, class U02,
           class U10, class U11, class U12,
           class U20, class U21, class U22>
  MatrixBlock_3x3( const std::initializer_list< std::initializer_list<U00> >& s00,
                   const std::initializer_list< std::initializer_list<U01> >& s01,
                   const std::initializer_list< std::initializer_list<U02> >& s02,

                   const std::initializer_list< std::initializer_list<U10> >& s10,
                   const std::initializer_list< std::initializer_list<U11> >& s11,
                   const std::initializer_list< std::initializer_list<U12> >& s12,

                   const std::initializer_list< std::initializer_list<U20> >& s20,
                   const std::initializer_list< std::initializer_list<U21> >& s21,
                   const std::initializer_list< std::initializer_list<U22> >& s22) :
    m00(s00), m01(s01), m02(s02),
    m10(s10), m11(s11), m12(s12),
    m20(s20), m21(s21), m22(s22) {}

  MatrixBlock_3x3( const typename Matrix00::size_type& s00, const typename Matrix01::size_type& s01, const typename Matrix02::size_type& s02,
                   const typename Matrix10::size_type& s10, const typename Matrix11::size_type& s11, const typename Matrix12::size_type& s12,
                   const typename Matrix20::size_type& s20, const typename Matrix21::size_type& s21, const typename Matrix22::size_type& s22) :
    m00(s00), m01(s01), m02(s02),
    m10(s10), m11(s11), m12(s12),
    m20(s20), m21(s21), m22(s22) {}

  // Generic constructor.
  template<class A00, class A01, class A02,
           class A10, class A11, class A12,
           class A20, class A21, class A22>
  MatrixBlock_3x3( const A00& a00, const A01& a01, const A02& a02,
                   const A10& a10, const A11& a11, const A12& a12,
                   const A20& a20, const A21& a21, const A22& a22 ) :
    m00(a00), m01(a01), m02(a02),
    m10(a10), m11(a11), m12(a12),
    m20(a20), m21(a21), m22(a22) {}


  // Generic constructor. For example, non-zero patterns used to construct sparse matricies
  template<class A00, class A01, class A02,
           class A10, class A11, class A12,
           class A20, class A21, class A22>
  MatrixBlock_3x3( const MatrixBlock_3x3<A00, A01, A02,
                                         A10, A11, A12,
                                         A20, A21, A22>& A ) :
    m00(A.m00), m01(A.m01), m02(A.m02),
    m10(A.m10), m11(A.m11), m12(A.m12),
    m20(A.m20), m21(A.m21), m22(A.m22) {}

  //Assignment operators
  MatrixBlock_3x3& operator=(const Real& v);

  MatrixBlock_3x3& operator=(const MatrixBlock_3x3& A);

  MatrixBlock_3x3& operator=( const DLA::Identity& I);

  template<class T>
  MatrixBlock_3x3& operator*=( const T& s );

  // Lazy expression functions
  template<class Vector0, class Vector1, class Vector2>
  void mulVec_value(const VectorBlock_3<Vector0, Vector1, Vector2>& x,
                    const Real sgn,
                    VectorBlock_3<Vector0, Vector1, Vector2>& b) const;

  template<class Vector0, class Vector1, class Vector2>
  void mulVec_plus(const VectorBlock_3<Vector0, Vector1, Vector2>& x,
                   const Real sgn,
                   VectorBlock_3<Vector0, Vector1, Vector2>& b) const;

  int m() const { return 3; }
  int n() const { return 3; }

  Matrix00 m00;
  Matrix01 m01;
  Matrix02 m02;

  Matrix10 m10;
  Matrix11 m11;
  Matrix12 m12;

  Matrix20 m20;
  Matrix21 m21;
  Matrix22 m22;
};

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
inline MatrixBlock_3x3<Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22>&
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22 >::
operator=(const Real& v)
{
  m00 = v; m01 = v; m02 = v;
  m10 = v; m11 = v; m12 = v;
  m20 = v; m21 = v; m22 = v;

  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
inline MatrixBlock_3x3<Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22>&
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22 >::
operator=(const MatrixBlock_3x3& A)
{
  m00 = A.m00;
  m01 = A.m01;
  m02 = A.m02;
  m10 = A.m10;
  m11 = A.m11;
  m12 = A.m12;
  m20 = A.m20;
  m21 = A.m21;
  m22 = A.m22;
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
inline MatrixBlock_3x3<Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22>&
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22 >::
operator=(const DLA::Identity& I)
{
  m00 = I; m01 = 0; m02 = 0;
  m10 = 0; m11 = I; m12 = 0;
  m20 = 0; m21 = 0; m22 = I;
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
template<class T>
inline MatrixBlock_3x3<Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22>&
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02, Matrix10, Matrix11, Matrix12, Matrix20, Matrix21, Matrix22 >::
operator*=(const T& s)
{
  m00 *= s; m01 *= s; m02 *= s;
  m10 *= s; m11 *= s; m12 *= s;
  m20 *= s; m21 *= s; m22 *= s;
  return *this;
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
template<class Vector0, class Vector1, class Vector2>
void
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02,
                 Matrix10, Matrix11, Matrix12,
                 Matrix20, Matrix21, Matrix22>::mulVec_value(const VectorBlock_3<Vector0, Vector1, Vector2>& x,
                                                             const Real sgn,
                                                             VectorBlock_3<Vector0, Vector1, Vector2>& b) const
{
  b = 0;
  mulVec_plus(x, sgn, b);
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02,
         class Matrix10, class Matrix11, class Matrix12,
         class Matrix20, class Matrix21, class Matrix22>
template<class Vector0, class Vector1, class Vector2>
void
MatrixBlock_3x3< Matrix00, Matrix01, Matrix02,
                 Matrix10, Matrix11, Matrix12,
                 Matrix20, Matrix21, Matrix22>::mulVec_plus(const VectorBlock_3<Vector0, Vector1, Vector2>& x,
                                                            const Real sgn,
                                                            VectorBlock_3<Vector0, Vector1, Vector2>& b) const
{
  b.v0 += sgn*(m00*x.v0 + m01*x.v1 + m02*x.v2);
  b.v1 += sgn*(m10*x.v0 + m11*x.v1 + m12*x.v2);
  b.v2 += sgn*(m20*x.v0 + m21*x.v1 + m22*x.v2);
}

} //namespace BLA
} //namespace numpack 

#endif //MATRIXBLOCK_2X2_H
