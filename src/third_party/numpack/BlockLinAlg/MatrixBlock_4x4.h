// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SRC_LINEARALGEBRA_BLOCKLINALG_MATRIXBLOCK_4X4_H_
#define SRC_LINEARALGEBRA_BLOCKLINALG_MATRIXBLOCK_4X4_H_

#include <initializer_list>

#include "BlockLinAlg_Mul.h"
#include "BlockLinAlg_Add.h"

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include "numpack/DenseLinAlg/tools/Identity.h"

//----------------------------------------------------------------------------//
// A class for generically storing 4x4 block matrices where the
// Matrices can be of mixed types
//
//----------------------------------------------------------------------------//

namespace numpack 
{

namespace BLA
{

template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
class MatrixBlock_4x4 : public BlockMatrixType< MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                                                                Matrix10, Matrix11, Matrix12, Matrix13,
                                                                Matrix20, Matrix21, Matrix22, Matrix23,
                                                                Matrix30, Matrix31, Matrix32, Matrix33> >,
                        public BlockLinAlgType< MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                                                                Matrix10, Matrix11, Matrix12, Matrix13,
                                                                Matrix20, Matrix21, Matrix22, Matrix23,
                                                                Matrix30, Matrix31, Matrix32, Matrix33> >
{
public:

#if 0 // not unit tested yet
  typedef MatrixBlock_4x4<typename Matrix00::size_type, typename Matrix01::size_type, typename Matrix02::size_type, typename Matrix03::size_type,
                          typename Matrix10::size_type, typename Matrix11::size_type, typename Matrix12::size_type, typename Matrix13::size_type,
                          typename Matrix20::size_type, typename Matrix21::size_type, typename Matrix22::size_type, typename Matrix23::size_type,
                          typename Matrix30::size_type, typename Matrix31::size_type, typename Matrix32::size_type, typename Matrix33::size_type>
          size_type;

  // Constructor: assemble individual initializer lists
  template<class U00, class U01, class U02, class U03,
           class U10, class U11, class U12, class U13,
           class U20, class U21, class U22, class U23,
           class U30, class U31, class U32, class U33>
  MatrixBlock_4x4( const std::initializer_list<U00>& s00, const std::initializer_list<U01>& s01,
                   const std::initializer_list<U02>& s02, const std::initializer_list<U03>& s03,
                   //
                   const std::initializer_list<U10>& s10, const std::initializer_list<U11>& s11,
                   const std::initializer_list<U12>& s12, const std::initializer_list<U13>& s13,
                   //
                   const std::initializer_list<U20>& s20, const std::initializer_list<U21>& s21,
                   const std::initializer_list<U22>& s22, const std::initializer_list<U23>& s23,
                   //
                   const std::initializer_list<U30>& s30, const std::initializer_list<U31>& s31,
                   const std::initializer_list<U32>& s32, const std::initializer_list<U33>& s33 ) :
    m00(s00), m01(s01), m02(s02), m03(s03),
    m10(s10), m11(s11), m12(s12), m13(s13),
    m20(s20), m21(s21), m22(s22), m23(s23),
    m30(s30), m31(s31), m32(s32), m33(s33) {}
#endif

  template<class U00, class U01, class U02, class U03,
           class U10, class U11, class U12, class U13,
           class U20, class U21, class U22, class U23,
           class U30, class U31, class U32, class U33>
  MatrixBlock_4x4( const std::initializer_list< std::initializer_list<U00> >& s00, const std::initializer_list< std::initializer_list<U01> >& s01,
                   const std::initializer_list< std::initializer_list<U02> >& s02, const std::initializer_list< std::initializer_list<U03> >& s03,
                   //
                   const std::initializer_list< std::initializer_list<U10> >& s10, const std::initializer_list< std::initializer_list<U11> >& s11,
                   const std::initializer_list< std::initializer_list<U12> >& s12, const std::initializer_list< std::initializer_list<U13> >& s13,
                   //
                   const std::initializer_list< std::initializer_list<U20> >& s20, const std::initializer_list< std::initializer_list<U21> >& s21,
                   const std::initializer_list< std::initializer_list<U22> >& s22, const std::initializer_list< std::initializer_list<U23> >& s23,
                   //
                   const std::initializer_list< std::initializer_list<U30> >& s30, const std::initializer_list< std::initializer_list<U31> >& s31,
                   const std::initializer_list< std::initializer_list<U32> >& s32, const std::initializer_list< std::initializer_list<U33> >& s33) :
    m00(s00), m01(s01), m02(s02), m03(s03),
    m10(s10), m11(s11), m12(s12), m13(s13),
    m20(s20), m21(s21), m22(s22), m23(s23),
    m30(s30), m31(s31), m32(s32), m33(s33) {}
#if 0 // not unit tested yet
  MatrixBlock_4x4( const typename Matrix00::size_type& s00, const typename Matrix01::size_type& s01,
                   const typename Matrix02::size_type& s02, const typename Matrix03::size_type& s03,
                   //
                   const typename Matrix10::size_type& s10, const typename Matrix11::size_type& s11,
                   const typename Matrix12::size_type& s12, const typename Matrix13::size_type& s13,
                   //
                   const typename Matrix20::size_type& s20, const typename Matrix21::size_type& s21,
                   const typename Matrix22::size_type& s22, const typename Matrix23::size_type& s23,
                   //
                   const typename Matrix30::size_type& s30, const typename Matrix31::size_type& s31,
                   const typename Matrix32::size_type& s32, const typename Matrix33::size_type& s33) :
    m00(s00), m01(s01), m02(s02), m03(s03),
    m10(s10), m11(s11), m12(s12), m13(s13),
    m20(s20), m21(s21), m22(s22), m23(s23),
    m30(s30), m31(s31), m32(s32), m33(s33) {}
#endif

  // TODO: unit test missing
  // Constructor: assemble generic individual blocks
  template<class A00, class A01, class A02, class A03,
           class A10, class A11, class A12, class A13,
           class A20, class A21, class A22, class A23,
           class A30, class A31, class A32, class A33>
  MatrixBlock_4x4( const A00& a00, const A01& a01, const A02& a02, const A03& a03,
                   const A10& a10, const A11& a11, const A12& a12, const A13& a13,
                   const A20& a20, const A21& a21, const A22& a22, const A23& a23,
                   const A30& a30, const A31& a31, const A32& a32, const A33& a33 ) :
    m00(a00), m01(a01), m02(a02), m03(a03),
    m10(a10), m11(a11), m12(a12), m13(a13),
    m20(a20), m21(a21), m22(a22), m23(a23),
    m30(a30), m31(a31), m32(a32), m33(a33) {}

  // Constructor: copy size-4x4 block matrix
  // For example, non-zero patterns used to construct block matricies
  template<class A00, class A01, class A02, class A03,
           class A10, class A11, class A12, class A13,
           class A20, class A21, class A22, class A23,
           class A30, class A31, class A32, class A33>
  MatrixBlock_4x4( const MatrixBlock_4x4<A00,A01,A02,A03,
                                         A10,A11,A12,A13,
                                         A20,A21,A22,A23,
                                         A30,A31,A32,A33>& A ) :
    m00(A.m00), m01(A.m01), m02(A.m02), m03(A.m03),
    m10(A.m10), m11(A.m11), m12(A.m12), m13(A.m13),
    m20(A.m20), m21(A.m21), m22(A.m22), m23(A.m23),
    m30(A.m30), m31(A.m31), m32(A.m32), m33(A.m33) {}

  // assignment operators
  MatrixBlock_4x4& operator=(const Real& v)
  {
    m00=v; m01=v; m02=v; m03=v;
    m10=v; m11=v; m12=v; m13=v;
    m20=v; m21=v; m22=v; m23=v;
    m30=v; m31=v; m32=v; m33=v;
    return *this;
  }

#if 0 // not working now because the copy assignment operator of A.m00 is implicitly deleted
  MatrixBlock_4x4& operator=(const MatrixBlock_4x4& A)
  {
    m00=A.m00; m01=A.m01; m02=A.m02; m03=A.m03;
    m10=A.m10; m11=A.m11; m12=A.m12; m13=A.m13;
    m20=A.m20; m21=A.m21; m22=A.m22; m23=A.m23;
    m30=A.m30; m31=A.m31; m32=A.m32; m33=A.m33;
    return *this;
  }
#endif

  MatrixBlock_4x4& operator=(const DLA::Identity& I)
  {
    m00=I; m01=0; m02=0; m03=0;
    m10=0; m11=I; m12=0; m13=0;
    m20=0; m21=0; m22=I; m23=0;
    m30=0; m31=0; m32=0; m33=I;
    return *this;
  }

#if 0 // not unit tested yet
  template<class Expr> MatrixBlock_4x4& operator=(const BlockLinAlgType<Expr>& expr)
  {
    const Expr& Tree = expr.cast();
    Tree.value(1,*this);
    return *this;
  }

  // binary assignment/accumulation
  inline void value(const Real sgn, MatrixBlock_4x4& A) const;
  inline void plus(const Real sgn, MatrixBlock_4x4& A) const;

  // Lazy evaluation expression/overloaded operators
  template<class Expr> inline MatrixBlock_4x4& operator+=(const BlockLinAlgType<Expr>& expr)
  {
    const Expr& Tree = expr.cast();
    Tree.plus(1,*this);
    return *this;
  }
  template<class Expr> inline MatrixBlock_4x4& operator-=(const BlockLinAlgType<Expr>& expr)
  {
    const Expr& Tree = expr.cast();
    Tree.plus(-1,*this);
    return *this;
  }
#if 0 // not sure if this is needed
  template<class Expr> MatrixBlock_4x4& operator*=( const BlockLinAlgType<Expr>& );
#endif

  template<class T>
  MatrixBlock_4x4& operator*=(const T& s)
  {
    m00*=s; m01*=s; m02*=s; m03*=s;
    m10*=s; m11*=s; m12*=s; m13*=s;
    m20*=s; m21*=s; m22*=s; m23*=s;
    m30*=s; m31*=s; m32*=s; m33*=s;
    return *this;
  }

#if 0 // working on it

  template<class Ta>
  inline void scale_row(const int i, const Ta& a);

  template<class Ta>
  inline void axpy_rows(const int i, const int j, const Ta& a);

  // Lazy expression functions
  template<class Vector0, class Vector1>
  void mulVec_value(const VectorBlock_2<Vector0, Vector1>& x,
                    const Real sgn,
                    VectorBlock_2<Vector0, Vector1>& b) const;

  template<class Vector0, class Vector1>
  void mulVec_plus(const VectorBlock_2<Vector0, Vector1>& x,
                   const Real sgn,
                   VectorBlock_2<Vector0, Vector1>& b) const;

  template<class M00, class M01, class M10, class M11>
  void mulMat_value(const MatrixBlock_4x4<M00, M01, M10, M11>& M,
                    const Real sgn,
                    MatrixBlock_4x4<M00, M01, M10, M11>& res) const;

  template<class M00, class M01, class M10, class M11>
  void mulMat_plus(const MatrixBlock_4x4<M00, M01, M10, M11>& M,
                   const Real sgn,
                   MatrixBlock_4x4<M00, M01, M10, M11>& res) const;
#endif
#endif

  // Lazy expression functions
  template<class Vector0, class Vector1, class Vector2, class Vector3>
  void mulVec_value(const VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& x,
                    const Real sgn,
                    VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& b) const;

  template<class Vector0, class Vector1, class Vector2, class Vector3>
  void mulVec_plus(const VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& x,
                   const Real sgn,
                   VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& b) const;

  // get dimension
  static int m() { return 4; }
  static int n() { return 4; }

  // block component data
  Matrix00 m00; Matrix01 m01; Matrix02 m02; Matrix03 m03;
  Matrix10 m10; Matrix11 m11; Matrix12 m12; Matrix13 m13;
  Matrix20 m20; Matrix21 m21; Matrix22 m22; Matrix23 m23;
  Matrix30 m30; Matrix31 m31; Matrix32 m32; Matrix33 m33;
};

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
template<class Vector0, class Vector1, class Vector2, class Vector3>
void
MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                Matrix10, Matrix11, Matrix12, Matrix13,
                Matrix20, Matrix21, Matrix22, Matrix23,
                Matrix30, Matrix31, Matrix32, Matrix33>::
mulVec_value(const VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& x,
             const Real sgn,
             VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& b) const
{
  b = 0;
  mulVec_plus(x, sgn, b);
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
template<class Vector0, class Vector1, class Vector2, class Vector3>
void
MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                Matrix10, Matrix11, Matrix12, Matrix13,
                Matrix20, Matrix21, Matrix22, Matrix23,
                Matrix30, Matrix31, Matrix32, Matrix33>::
mulVec_plus(const VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& x,
            const Real sgn,
            VectorBlock_4<Vector0, Vector1, Vector2, Vector3>& b) const
{
  b.v0 += sgn*(m00*x.v0 + m01*x.v1 + m02*x.v2 + m03*x.v3);
  b.v1 += sgn*(m10*x.v0 + m11*x.v1 + m12*x.v2 + m13*x.v3);
  b.v2 += sgn*(m20*x.v0 + m21*x.v1 + m22*x.v2 + m23*x.v3);
  b.v3 += sgn*(m30*x.v0 + m31*x.v1 + m32*x.v2 + m33*x.v3);
}

#if 0 // not unit tested yet
//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
void
MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                Matrix10, Matrix11, Matrix12, Matrix13,
                Matrix20, Matrix21, Matrix22, Matrix23,
                Matrix30, Matrix31, Matrix32, Matrix33>::
value(const Real sgn, MatrixBlock_4x4& res) const
{
  res.m00=sgn*m00; res.m01=sgn*m01; res.m02=sgn*m02; res.m03=sgn*m03;
  res.m10=sgn*m10; res.m11=sgn*m11; res.m12=sgn*m12; res.m13=sgn*m13;
  res.m20=sgn*m20; res.m21=sgn*m21; res.m22=sgn*m22; res.m23=sgn*m23;
  res.m30=sgn*m30; res.m31=sgn*m31; res.m32=sgn*m32; res.m33=sgn*m33;
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
void
MatrixBlock_4x4<Matrix00, Matrix01, Matrix02, Matrix03,
                Matrix10, Matrix11, Matrix12, Matrix13,
                Matrix20, Matrix21, Matrix22, Matrix23,
                Matrix30, Matrix31, Matrix32, Matrix33>::
plus(const Real sgn, MatrixBlock_4x4& res) const
{
  res.m00+=sgn*m00; res.m01+=sgn*m01; res.m02+=sgn*m02; res.m03+=sgn*m03;
  res.m10+=sgn*m10; res.m11+=sgn*m11; res.m12+=sgn*m12; res.m13+=sgn*m13;
  res.m20+=sgn*m20; res.m21+=sgn*m21; res.m22+=sgn*m22; res.m23+=sgn*m23;
  res.m30+=sgn*m30; res.m31+=sgn*m31; res.m32+=sgn*m32; res.m33+=sgn*m33;
}

#if 0 // working on it
//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
template<class M00, class M01, class M10, class M11>
void
MatrixBlock_4x4< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulMat_value(const MatrixBlock_4x4<M00, M01, M10, M11>& M,
                                                    const Real sgn,
                                                    MatrixBlock_4x4<M00, M01, M10, M11>& res) const
{
  res = 0;
  mulMat_plus(M, sgn, res);
}

//=============================================================================
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
template<class M00, class M01, class M10, class M11>
void
MatrixBlock_4x4< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::mulMat_plus(const MatrixBlock_4x4<M00, M01, M10, M11>& M,
                                                   const Real sgn,
                                                   MatrixBlock_4x4<M00, M01, M10, M11>& res) const
{
  res.m00 += sgn*(m00*M.m00 + m01*M.m10);
  res.m01 += sgn*(m00*M.m01 + m01*M.m11);
  res.m10 += sgn*(m10*M.m00 + m11*M.m10);
  res.m11 += sgn*(m10*M.m01 + m11*M.m11);
}

//---------------------------------------------------------------------------//
template<class Matrix00, class Matrix01, class Matrix02, class Matrix03,
         class Matrix10, class Matrix11, class Matrix12, class Matrix13,
         class Matrix20, class Matrix21, class Matrix22, class Matrix23,
         class Matrix30, class Matrix31, class Matrix32, class Matrix33>
template<class Ta>
inline void
MatrixBlock_4x4< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::scale_row(const int i, const Ta& a)
{
  if (i == 0)
  {
    Matrix00 tmp00 = a*m00;
    Matrix01 tmp01 = a*m01;
    m00 = tmp00;
    m01 = tmp01;
  }
  else if (i == 1)
  {
    Matrix10 tmp10 = a*m10;
    Matrix11 tmp11 = a*m11;
    m10 = tmp10;
    m11 = tmp11;
  }
  else
    SANS_DEVELOPER_EXCEPTION("MatrixBlock_4x4::scale_row - row index i (=%d) has to be 0 or 1.", i);
}

//Performs y = a*x + y where x = row(i) and y = row(j)
template<class Matrix00, class Matrix01,
         class Matrix10, class Matrix11>
template<class Ta>
inline void
MatrixBlock_4x4< Matrix00, Matrix01,
                 Matrix10, Matrix11 >::axpy_rows(const int i, const int j, const Ta& a)
{
  if (i == 0)
  {
    if (j == 0)
    {
      m00 += a*m00;
      m01 += a*m01;
    }
    else if (j == 1)
    {
      m10 += a*m00;
      m11 += a*m01;
    }
    else
      SANS_DEVELOPER_EXCEPTION("MatrixBlock_4x4::axpy_rows - row index j (=%d) has to be 0 or 1.", j);
  }
  else if (i == 1)
  {
    if (j == 0)
    {
      m00 += a*m10;
      m01 += a*m11;
    }
    else if (j == 1)
    {
      m10 += a*m10;
      m11 += a*m11;
    }
    else
      SANS_DEVELOPER_EXCEPTION("MatrixBlock_4x4::axpy_rows - row index j (=%d) has to be 0 or 1.", j);
  }
  else
    SANS_DEVELOPER_EXCEPTION("MatrixBlock_4x4::axpy_rows - row index i (=%d) has to be 0 or 1.", i);
}
#endif

#endif
} //namespace BLA
} //namespace numpack 

#endif /* SRC_LINEARALGEBRA_BLOCKLINALG_MATRIXBLOCK_4X4_H_ */
