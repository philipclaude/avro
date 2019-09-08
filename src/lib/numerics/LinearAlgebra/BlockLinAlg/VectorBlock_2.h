// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORBLOCK_2_H
#define VECTORBLOCK_2_H

#include <initializer_list>
#include <iostream>

#include "tools/SANSnumerics.h"
#include "tools/SANSException.h"

#include "BlockLinAlg_Type.h"

//----------------------------------------------------------------------------//
// A class for generically storing 2 block vector where the
// vectors can be of mixed types
//
//----------------------------------------------------------------------------//

namespace SANS
{

namespace BLA
{

template<class Vector0,
         class Vector1>
class VectorBlock_2 : public BlockVectorType< VectorBlock_2<Vector0, Vector1> >,
                      public BlockLinAlgType< VectorBlock_2<Vector0, Vector1> >
{
public:

  template<class U0,
           class U1>
  VectorBlock_2( const std::initializer_list<U0>& s0,
                 const std::initializer_list<U1>& s1 ) :
    v0(s0),
    v1(s1) {}

  // Generic constructor.
  template<class V0,
           class V1>
  VectorBlock_2( const V0& v0,
                 const V1& v1 ) :
    v0(v0),
    v1(v1) {}


  // Generic constructor.
  template<class V0,
           class V1>
  VectorBlock_2( const VectorBlock_2<V0,
                                     V1>& V ) :
    v0(V.v0),
    v1(V.v1) {}

  VectorBlock_2& operator=(const Real& v) { v0 = v; v1 = v; return *this; }
  VectorBlock_2& operator=(const VectorBlock_2& v) { v0 = v.v0; v1 = v.v1; return *this; }

  template<class T>
  VectorBlock_2& operator*=( const T& s ) { v0 *= s; v1 *= s; return *this; }

  //Lazy evaluation functions
  template<class V0, class V1>
  void value(const Real sgn, VectorBlock_2<V0, V1>& b) const;
  template<class V0, class V1>
  void plus(const Real sgn, VectorBlock_2<V0, V1>& b) const;

  template< class Expr > inline VectorBlock_2& operator=( const BlockLinAlgType<Expr>& r );
  template< class Expr > inline VectorBlock_2& operator+=( const BlockLinAlgType<Expr>& r );
  template< class Expr > inline VectorBlock_2& operator-=( const BlockLinAlgType<Expr>& r );

  template<class Ta>
  inline void scale_row0(const Ta& a);
  template<class Ta>
  inline void scale_row1(const Ta& a);

  template<class Ta>
  inline void axpy_rows10(const Ta& a);
  template<class Ta>
  inline void axpy_rows01(const Ta& a);

  int m() const { return 2; }
  int n() const { return 1; }

  Vector0 v0;
  Vector1 v1;
};

//---------------------------------------------------------------------------//
template<class Vector0,
         class Vector1>
template< class V0, class V1 >
inline void
VectorBlock_2< Vector0, Vector1 >::value(const Real sgn, VectorBlock_2<V0, V1>& b) const
{
  b.v0 = sgn*v0;
  b.v1 = sgn*v1;
}

//---------------------------------------------------------------------------//
template<class Vector0,
         class Vector1>
template< class V0, class V1 >
inline void
VectorBlock_2< Vector0, Vector1 >::plus(const Real sgn, VectorBlock_2<V0, V1>& b) const
{
  b.v0 += sgn*v0;
  b.v1 += sgn*v1;
}


//---------------------------------------------------------------------------//
//Lazy assignment operators
template<class Vector0,
         class Vector1>
template< class Expr >
inline VectorBlock_2< Vector0, Vector1 >&
VectorBlock_2< Vector0, Vector1 >::operator=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.value(1, *this);

  return *this;
}

template<class Vector0,
         class Vector1>
template< class Expr >
inline VectorBlock_2< Vector0, Vector1 >&
VectorBlock_2< Vector0, Vector1 >::operator+=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.plus(1, *this);

  return *this;
}

template<class Vector0,
         class Vector1>
template< class Expr >
inline VectorBlock_2< Vector0, Vector1 >&
VectorBlock_2< Vector0, Vector1 >::operator-=( const BlockLinAlgType<Expr>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m() == Tree.m(), "with m() = %d and Tree.m() = %d", m(), Tree.m() );

  Tree.plus(-1, *this);

  return *this;
}

template<class Vector0,
         class Vector1>
template<class Ta>
inline void
VectorBlock_2< Vector0, Vector1 >::scale_row0(const Ta& a)
{
  Vector0 tmp(v0.size());
  tmp = a*v0;
  v0 = tmp;
}


template<class Vector0,
         class Vector1>
template<class Ta>
inline void
VectorBlock_2< Vector0, Vector1 >::scale_row1(const Ta& a)
{
  Vector1 tmp(v1.size());
  tmp = a*v1;
  v1 = tmp;
}


//Performs y = a*x + y where x = row(0) and y = row(1)
template<class Vector0,
         class Vector1>
template<class Ta>
inline void
VectorBlock_2< Vector0, Vector1 >::axpy_rows01(const Ta& a)
{
  v1 += a*v0;
}

//Performs y = a*x + y where x = row(1) and y = row(0)
template<class Vector0,
         class Vector1>
template<class Ta>
inline void
VectorBlock_2< Vector0, Vector1 >::axpy_rows10(const Ta& a)
{
  v0 += a*v1;
}

// I/O
template<class Vector0, class Vector1>
std::ostream&
operator<<( std::ostream& out, const VectorBlock_2<Vector0,Vector1>& v )
{
  out << "VectorBlock_2.v0 = " << std::endl;
  out << v.v0 << std::endl;
  out << std::endl;

  out << "VectorBlock_2.v1 = " << std::endl;
  out << v.v1 << std::endl;
  out << std::endl;

  return out;
}

} //namespace BLA
} //namespace SANS

#endif //VECTORBLOCK_2_H
