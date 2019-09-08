// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORBLOCK_3_H
#define VECTORBLOCK_3_H

#include <initializer_list>

#include "tools/SANSnumerics.h"
#include "BlockLinAlg_Type.h"

//----------------------------------------------------------------------------//
// A class for generically storing 3 block vector where the
// vectors can be of mixed types
//
//----------------------------------------------------------------------------//

namespace SANS
{

namespace BLA
{

template<class Vector0,
         class Vector1,
         class Vector2>
class VectorBlock_3 : public BlockVectorType< VectorBlock_3<Vector0, Vector1, Vector2> >,
                      public BlockLinAlgType< VectorBlock_3<Vector0, Vector1, Vector2> >
{
public:

  template<class U0, class U1, class U2>
  VectorBlock_3( const std::initializer_list<U0>& s0,
                 const std::initializer_list<U1>& s1,
                 const std::initializer_list<U2>& s2) :
    v0(s0),
    v1(s1),
    v2(s2) {}

  // Generic constructor.
  template<class V0, class V1, class V2>
  VectorBlock_3( const V0& v0,
                 const V1& v1,
                 const V2& v2) :
    v0(v0),
    v1(v1),
    v2(v2) {}


  // Generic constructor.
  template<class V0, class V1, class V2>
  VectorBlock_3( const VectorBlock_3<V0,
                                     V1,
                                     V2>& V ) :
    v0(V.v0),
    v1(V.v1),
    v2(V.v2) {}

  VectorBlock_3& operator=(const Real& v) { v0 = v; v1 = v; v2 = v; return *this; }
  VectorBlock_3& operator=(const VectorBlock_3& v) { v0 = v.v0; v1 = v.v1; v2 = v.v2; return *this; }

  template<class T>
  VectorBlock_3& operator*=( const T& s ) { v0 *= s; v1 *= s; v2 *= s; return *this; }

  //Lazy evaluation functions
  template<class V0, class V1, class V2>
  void value(const Real sgn, VectorBlock_3<V0, V1, V2>& b) const;
  template<class V0, class V1, class V2>
  void plus(const Real sgn, VectorBlock_3<V0, V1, V2>& b) const;

  template< class Expr > inline VectorBlock_3& operator=( const BlockLinAlgType<Expr>& r );
  template< class Expr > inline VectorBlock_3& operator+=( const BlockLinAlgType<Expr>& r );
  template< class Expr > inline VectorBlock_3& operator-=( const BlockLinAlgType<Expr>& r );

  int m() const { return 3; }
  int n() const { return 1; }

  Vector0 v0;
  Vector1 v1;
  Vector2 v2;
};

//---------------------------------------------------------------------------//
template<class Vector0, class Vector1, class Vector2>
template< class V0, class V1, class V2 >
inline void
VectorBlock_3< Vector0, Vector1, Vector2 >::value(const Real sgn, VectorBlock_3<V0, V1, V2>& b) const
{
  b.v0 = sgn*v0;
  b.v1 = sgn*v1;
  b.v2 = sgn*v2;
}

//---------------------------------------------------------------------------//
template<class Vector0, class Vector1, class Vector2>
template< class V0, class V1, class V2 >
inline void
VectorBlock_3< Vector0, Vector1, Vector2 >::plus(const Real sgn, VectorBlock_3<V0, V1, V2>& b) const
{
  b.v0 += sgn*v0;
  b.v1 += sgn*v1;
  b.v2 += sgn*v2;
}


//---------------------------------------------------------------------------//
//Lazy assignment operators
template<class Vector0, class Vector1, class Vector2>
template< class Expr >
inline VectorBlock_3< Vector0, Vector1, Vector2 >&
VectorBlock_3< Vector0, Vector1, Vector2 >::operator=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.value(1, *this);

  return *this;
}

template<class Vector0, class Vector1, class Vector2>
template< class Expr >
inline VectorBlock_3< Vector0, Vector1, Vector2 >&
VectorBlock_3< Vector0, Vector1, Vector2 >::operator+=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.plus(1, *this);

  return *this;
}

template<class Vector0, class Vector1, class Vector2>
template< class Expr >
inline VectorBlock_3< Vector0, Vector1, Vector2 >&
VectorBlock_3< Vector0, Vector1, Vector2 >::operator-=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.plus(-1, *this);

  return *this;
}

} //namespace BLA
} //namespace SANS

#endif //VECTORBLOCK_3_H
