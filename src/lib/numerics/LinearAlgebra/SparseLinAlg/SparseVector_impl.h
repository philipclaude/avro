// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//No include block because this file needs to be included twice

#ifdef TV_SPEC


namespace SANS
{
namespace SLA
{

//=============================================================================
//Implementation of the SparseVector member functions
template<class TV>
SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator=(const Real v)
{
  for (int i = 0; i < m_; i++)
    (*this)[i] = v;

  return *this;
}

template<class TV>
SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator=(const SparseVector& v)
{
  SANS_ASSERT_MSG( m_ == v.m(), "with m_ = %d and v.m() = %d", m_, v.m() );
  for (int i = 0; i < m_; i++)
    (*this)[i] = v[i];

  return *this;
}

template<class TV>
SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator/=(const Real v)
{
  for (int i = 0; i < m_; i++)
    (*this)[i] /= v;

  return *this;
}

template<class TV>
SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator*=(const Real v)
{
  for (int i = 0; i < m_; i++)
    (*this)[i] *= v;

  return *this;
}

//Lazy evaluation functions
template<class TV>
void
SparseVector< TV_SPEC >::value(const Real sgn, SparseVector& b) const
{
  SANS_ASSERT( m_ == b.m() );

  for (int i = 0; i < m_; i++)
    b[i] = sgn*(*this)[i];
}

template<class TV>
void
SparseVector< TV_SPEC >::plus(const Real sgn, SparseVector& b) const
{
  SANS_ASSERT( m_ == b.m() );

  for (int i = 0; i < m_; i++)
    b[i] += sgn*(*this)[i];
}


//Lazy assignment operators
template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator=( const SparseLinAlgType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );

  Tree.value(1, *this);

  return *this;
}

template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator+=( const SparseLinAlgType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT_MSG( m_ == Tree.m(), "with m_ = %d and Tree.m() = %d", m_, Tree.m() );

  Tree.plus(1, *this);

  return *this;
}

template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator-=( const SparseLinAlgType<Expr, true>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );

  Tree.plus(-1, *this);

  return *this;
}

template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator=( const SparseLinAlgType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ could change in the loop
  const int tmp = m_;

  for ( int i = 0; i < tmp; i++)
    (*this)[i] = Tree[i];

  return *this;
}

template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator+=( const SparseLinAlgType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ could change in the loop
  const int tmp = m_;

  for ( int i = 0; i < tmp; i++)
    (*this)[i] += Tree[i];

  return *this;
}

template<class TV>
template< class Expr >
inline SparseVector< TV_SPEC >&
SparseVector< TV_SPEC >::operator-=( const SparseLinAlgType<Expr, false>& r )
{
  const Expr& Tree = r.cast();

  SANS_ASSERT( m_ == Tree.m() );

  //Saving to a local temporary improves optimization as the compiler does not
  //need to assume that m_ could change in the loop
  const int tmp = m_;

  for ( int i = 0; i < tmp; i++)
    (*this)[i] -= Tree[i];

  return *this;
}

} //namespace SLA
} //namespace SANS

#endif
