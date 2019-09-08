// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SRC_LINEARALGEBRA_BLOCKLINALG_VECTORBLOCK_4_H_
#define SRC_LINEARALGEBRA_BLOCKLINALG_VECTORBLOCK_4_H_

#include <initializer_list>
#include <iostream>

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include "BlockLinAlg_Type.h"

//----------------------------------------------------------------------------//
// A class for generically storing size-4 block vector where the
// vectors can be of mixed types (e.g. dense or sparse)
//
//----------------------------------------------------------------------------//

namespace SANS
{

namespace BLA
{

template<class Vector0, class Vector1, class Vector2, class Vector3>
class VectorBlock_4 : public BlockVectorType< VectorBlock_4<Vector0,Vector1,Vector2,Vector3> >,
                      public BlockLinAlgType< VectorBlock_4<Vector0,Vector1,Vector2,Vector3> >
{
public:
  // Constructor: assemble individual initializer lists
  template<class U0, class U1, class U2, class U3>
  VectorBlock_4( const std::initializer_list<U0>& s0,
                 const std::initializer_list<U1>& s1,
                 const std::initializer_list<U2>& s2,
                 const std::initializer_list<U3>& s3 ) :
    v0(s0), v1(s1), v2(s2), v3(s3) {}

  // Constructor: assemble generic individual vector block components
  template<class V0, class V1, class V2, class V3>
  VectorBlock_4(const V0& v0, const V1& v1, const V2& v2, const V3& v3) :
    v0(v0), v1(v1), v2(v2), v3(v3) {}

  // Move constructor: assemble generic individual vector block components
  template<class V0, class V1, class V2, class V3>
  VectorBlock_4(const V0&& v0, const V1&& v1, const V2&& v2, const V3&& v3) :
    v0(v0), v1(v1), v2(v2), v3(v3) {}

  // Constructor: copy size-4 block vector
  template<class V0, class V1, class V2, class V3>
  explicit VectorBlock_4(const VectorBlock_4<V0,V1,V2,V3>& V) :
    v0(V.v0), v1(V.v1), v2(V.v2), v3(V.v3) {}

  // Copy constructor
  VectorBlock_4(const VectorBlock_4& V) :
    v0(V.v0), v1(V.v1), v2(V.v2), v3(V.v3) {}

  // assignment operators
  VectorBlock_4& operator=(const Real& v) { v0 = v; v1 = v; v2 = v; v3 = v; return *this; }
  VectorBlock_4& operator=(const VectorBlock_4& v) { v0 = v.v0; v1 = v.v1; v2 = v.v2; v3 = v.v3; return *this; }
  template<class Expr> inline VectorBlock_4& operator=(const BlockLinAlgType<Expr>& r);

  // Lazy evaluation/assignment functions
  //
  // assign its value to b up to sign sgn
  template<class V0, class V1, class V2, class V3>
  void value(const Real sgn, VectorBlock_4<V0,V1,V2,V3>& b) const;
  // add its value to b up to sign sgn
  template<class V0, class V1, class V2, class V3>
  void plus(const Real sgn, VectorBlock_4<V0,V1,V2,V3>& b) const;

  template<class Expr> inline VectorBlock_4& operator+=(const BlockLinAlgType<Expr>& r);
  template<class Expr> inline VectorBlock_4& operator-=(const BlockLinAlgType<Expr>& r);

  template<class T>
  inline VectorBlock_4& operator*=(const T& s) { v0 *= s; v1 *= s; v2 *= s; v3 *= s; return *this; }

//  // scale a block row by a
//  template<class Ta>
//  inline void scale_row(const int i, const Ta& a);
//
//  template<class Ta>
//  inline void axpy_rows(const int i, const int j, const Ta& a);

  // get dimension
  static int m() { return 4; }
  static int n() { return 1; }

  Vector0 v0;
  Vector1 v1;
  Vector2 v2;
  Vector3 v3;
};

//---------------------------------------------------------------------------//
//Lazy assignment operators
//
template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class Expr>
inline VectorBlock_4<Vector0,Vector1,Vector2,Vector3>&
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
operator=(const BlockLinAlgType<Expr>& r)
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG(m()==Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m());

  Tree.value(+1, *this);

  return *this;
}

//---------------------------------------------------------------------------//
// Lazy evaluation/assignment functions
//
template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class V0, class V1, class V2, class V3>
inline void
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
value(const Real sgn, VectorBlock_4<V0,V1,V2,V3>& b) const
{
  b.v0 = sgn*v0;
  b.v1 = sgn*v1;
  b.v2 = sgn*v2;
  b.v3 = sgn*v3;
}

template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class V0, class V1, class V2, class V3>
inline void
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
plus(const Real sgn, VectorBlock_4<V0,V1,V2,V3>& b) const
{
  b.v0 += sgn*v0;
  b.v1 += sgn*v1;
  b.v2 += sgn*v2;
  b.v3 += sgn*v3;
}

template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class Expr>
inline VectorBlock_4<Vector0,Vector1,Vector2,Vector3>&
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
operator+=(const BlockLinAlgType<Expr>& r)
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG(m()==Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m());

  Tree.plus(+1, *this);

  return *this;
}

template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class Expr>
inline VectorBlock_4<Vector0,Vector1,Vector2,Vector3>&
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
operator-=(const BlockLinAlgType<Expr>& r)
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG(m()==Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.plus(-1, *this);

  return *this;
}

#if 0 // TODO: not needed for now
// Scale a block row by
template<class Vector0, class Vector1, class Vector2, class Vector3>
template<class Ta>
inline void
VectorBlock_4<Vector0,Vector1,Vector2,Vector3>::
scale_row(const int i, const Ta& a)
{
  if (i==0)
  {
    v0 *= a;
  }
  else if (i==1)
  {
    v1 *= a;
  }
  else if (i==2)
  {
    v2 *= a;
  }
  else if (i==3)
  {
    v3 *= a;
  }
  else
    SANS_DEVELOPER_EXCEPTION("VectorBlock_4<V0,V1,V2,V3>::scale_row - row index i (=%d) has to be in {0,1,2,3}", i);
}
#endif

// I/O
template<class Vector0, class Vector1, class Vector2, class Vector3>
std::ostream&
operator<<( std::ostream& out, const VectorBlock_4<Vector0,Vector1,Vector2,Vector3>& v )
{
  out << "VectorBlock_4.v0 = " << std::endl;
  out << v.v0 << std::endl;
  out << std::endl;

  out << "VectorBlock_4.v1 = " << std::endl;
  out << v.v1 << std::endl;
  out << std::endl;

  out << "VectorBlock_4.v2 = " << std::endl;
  out << v.v2 << std::endl;
  out << std::endl;

  out << "VectorBlock_4.v3 = " << std::endl;
  out << v.v3 << std::endl;
  out << std::endl;

  return out;
}

} //namespace BLA
} //namespace SANS

#endif /* SRC_LINEARALGEBRA_BLOCKLINALG_VECTORBLOCK_4_H_ */
