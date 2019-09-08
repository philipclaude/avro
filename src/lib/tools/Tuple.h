// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef TUPLE_H
#define TUPLE_H

#include <type_traits>
#include "abs.h"

namespace SANS
{

// This is a nice discussion on variadic templates
// http://eli.thegreenplace.net/2014/variadic-templates-in-c/

// A unique identifier used to indicate that a templated class is indeed a tuple
// The level provides a means of creating tuples of tuples of tuples
template<int Level = 0>
class TupleClass;

//=============================================================================
// This simplifies the declaration of a Tuple type, i.e.
//
// typedef Type0 Arg0;
// typedef Type1 Arg1;
// typedef Type2 Arg2;
// typedef typename MakeTuple<Tuple,Arg0,Arg1,Arg2>::type Args;
//
template<template<class,class,class> class Tuple, class... Args>
struct MakeTuple;

//Termination of the recursive loop: base case of variadic templates
template<template<class,class,class> class Tuple, class CompleteTuple>
struct MakeTuple<Tuple, CompleteTuple>
{
  typedef CompleteTuple type;
};

//Recursively build up the 0-Level Tuple
template<template<class,class,class> class Tuple, class L, class R, class... Args>
struct MakeTuple<Tuple, L, R, Args...>
{
  typedef typename MakeTuple< Tuple, Tuple<L, R, TupleClass<0> >, Args...>::type type;
};

//=============================================================================
// This simplifies the declaration of a Tuple where all arguments are the same type, i.e.
//
// const int N = 3;
// typedef typename MakeTupleN<Tuple,Type0,N>::type Args;
//
namespace detail
{

//Recursively build up the 0-Level Tuple
template<template<class,class,class> class Tuple, class TupleType, class Type, int N>
struct MakeTupleN_impl
{
  typedef typename MakeTupleN_impl< Tuple, Tuple<TupleType, Type, TupleClass<0> >, Type, N-1 >::type type;
};

//Termination of the recursive loop
template<template<class,class,class> class Tuple, class CompleteTuple, class Type>
struct MakeTupleN_impl<Tuple, CompleteTuple, Type, 0>
{
  typedef CompleteTuple type;
};

}

//Build up the Tuple
template<template<class,class,class> class Tuple, class Type, int N>
struct MakeTupleN
{
  static_assert( N >= 2, "The tuple must have at least 2 entries." );
  typedef typename detail::MakeTupleN_impl< Tuple, Tuple<Type, Type, TupleClass<0> >, Type, N-2 >::type type;
};

//=============================================================================
//A recursive template function to compute the number of arguments in a Tuple
template<class Type, int Level>
struct TupleSize
{
  static const int size = 1;
};

template< template<class,class,class> class Tuple, class L, class R, int Level>
struct TupleSize< Tuple<L, R, TupleClass<Level> >, Level >
{
  static const int size = 1 + TupleSize<L, Level>::size;
};

namespace detail
{
//=============================================================================
// This class locates the type associated with entry k in a Tuple.
// This is used in association with the get<k> function to find the appropriate return type
template<int k, int m, class Type, int Level = 0>
struct TupleType_impl;

template<int k, int m, template<class,class,class> class Tuple, class L, class R, int Level>
struct TupleType_impl<k, m, Tuple<L, R, TupleClass<Level> >, Level >
{
  typedef typename TupleType_impl<k,m-1,L,Level>::type type;
};

template<int k, template<class,class,class> class Tuple, class L, class R, int Level>
struct TupleType_impl<k, k, Tuple<L, R, TupleClass<Level> >, Level >
{
  typedef R type;
};

template<class L, int Level>
struct TupleType_impl<0,0,L,Level>
{
  typedef L type;
};
}

template<int k, class Type, int Level = 0>
struct TupleType
{
  static_assert( k == 0, "k must be zero if Type is not a Tuple");
  typedef Type type;
};

template<int k, template<class,class,class> class Tuple, class L, class R, int Level >
struct TupleType<k, Tuple<L, R, TupleClass<Level> >, Level >
{
  static const int size = TupleSize< Tuple<L,R,TupleClass<Level>>, Level >::size-1;
  typedef typename detail::TupleType_impl<k, size, Tuple<L,R,TupleClass<Level> >, Level >::type type;
};


//=============================================================================
// A class to allow for negative indexes extracting elements in re verse order
namespace detail
{

template<bool, int k, class Tuple, int Level>
struct TupleIndex;

template< int k, template<class,class,class> class Tuple, class L, class R, int Level>
struct TupleIndex<false, k, Tuple<L, R, TupleClass<Level>>, Level >
{
  static const int index = k;
};

template< int k, template<class,class,class> class Tuple, class L, class R, int Level>
struct TupleIndex<true, k, Tuple<L, R, TupleClass<Level>>, Level >
{
  static const int index = TupleSize< Tuple<L, R, TupleClass<Level>>, Level >::size + k;
};

}


//=============================================================================
// A simple interface to get the entry at index k in an instance of Tuple<L,R,TupleClass<> >, e.g.
//
// const Type0& t0 = get<0>(tuple);
// const Type1& t1 = get<1>(tuple);
// const Type2& t2 = get<2>(tuple);
//
namespace detail
{

template <int k, int m, int Level>
struct TupleGet_impl
{
  template < template<class,class,class> class Tuple, class L, class R>
  static const typename TupleType<k, Tuple<L,R,TupleClass<Level>>, Level >::type&
  get(const Tuple<L, R, TupleClass<Level>>& tuple)
  {
    return TupleGet_impl<k,m-1,Level>::get(tuple.left());
  }

  template < template<class,class,class> class Tuple, class L, class R>
  static typename TupleType<k, Tuple<L,R,TupleClass<Level>>, Level >::type&
  set(Tuple<L, R, TupleClass<Level>>& tuple)
  {
    return TupleGet_impl<k,m-1,Level>::set(tuple.left());
  }
};

template <int k, int Level>
struct TupleGet_impl<k, k, Level>
{
  template < template<class,class,class> class Tuple, class L, class R>
  static const typename TupleType<k,Tuple<L,R,TupleClass<Level>>, Level >::type&
  get(const Tuple<L, R, TupleClass<Level>>& tuple)
  {
    return tuple.right();
  }

  template < template<class,class,class> class Tuple, class L, class R>
  static typename TupleType<k,Tuple<L,R,TupleClass<Level>>, Level >::type&
  set(Tuple<L, R, TupleClass<Level>>& tuple)
  {
    return tuple.right();
  }
};

template <int Level>
struct TupleGet_impl<0, 0, Level>
{
  template <class L>
  static const L&
  get(const L& left)
  {
    return left;
  }

  template <class L>
  static L&
  set(L& left)
  {
    return left;
  }
};

}

template <int k, template<class,class,class> class Tuple, class L, class R, int Level>
const typename TupleType< detail::TupleIndex< (k < 0), k, Tuple<L,R,TupleClass<Level>>, Level >::index, Tuple<L,R,TupleClass<Level>>, Level >::type&
get(const Tuple<L, R, TupleClass<Level>>& tuple)
{
  static const int index = detail::TupleIndex< (k < 0), k, Tuple<L,R,TupleClass<Level>>, Level >::index;
  static const int size = TupleSize< Tuple<L, R, TupleClass<Level>>, Level >::size;
  static_assert(  size >  k, "k must be less than the number of elements" );
  static_assert( -size <= k, "negative k must be greater or equal to the negative size" );
  return detail::TupleGet_impl<index,size-1,Level>::get(tuple);
}

template <int k, template<class,class,class> class Tuple, class L, class R, int Level>
typename TupleType< detail::TupleIndex< (k < 0), k, Tuple<L,R,TupleClass<Level>>, Level >::index, Tuple<L,R,TupleClass<Level>>, Level >::type&
set(Tuple<L, R, TupleClass<Level>>& tuple)
{
  static const int index = detail::TupleIndex< (k < 0), k, Tuple<L,R,TupleClass<Level>>, Level >::index;
  static const int size = TupleSize< Tuple<L, R, TupleClass<Level>>, Level >::size;
  static_assert(  size >  k, "positive k must be less than the number of elements" );
  static_assert( -size <= k, "negative k must be greater or equal to the negative size" );
  return detail::TupleGet_impl<index,size-1,Level>::set(tuple);
}


//=============================================================================
// A simple interface to retrieve a specific Type in an instance of Tuple<L,R,TupleClass<Level> >, e.g.
//
// const double&       t0 = get<double>(tuple);
// const float&        t1 = get<float>(tuple);
// const ElementClass& t2 = get<ElementClass>(tuple);
//
// The first occurrence of the "Type" starting from the right that occurs in the tuple will be returned.
//
// Note that L and R are also allowed to be Tuple of levels lower than Level, e.g.
//
// Tuple< Tuple<L,R,TupleClass<LevelL> > , Tuple<L,R,TupleClass<LevelR> >, TupleClass<Level> >
//
// If no "Type" is found inside the current Tuple, then get() will try to find "Type" in the leftmost Tuple component of the original Tuple.
//

namespace detail
{

template <class Type>
struct TupleTypeGet_impl
{
  template < template<class,class,class> class Tuple, class L, class R, int Level>
  static const Type&
  get(const Tuple<L, R, TupleClass<Level> >& tuple)
  {
    return get(tuple.left());
  }

  template < template<class,class,class> class Tuple, class L, int Level>
  static const Type&
  get(const Tuple<L, Type, TupleClass<Level> >& tuple)
  {
    return tuple.right();
  }

  static const Type&
  get(const Type& type)
  {
    return type;
  }
};

}

template <class Type, template<class,class,class> class Tuple, class L, class R, int Level>
const Type&
get(const Tuple<L, R, TupleClass<Level> >& tuple)
{
  return detail::TupleTypeGet_impl<Type>::get(tuple);
}

}

#endif //TUPLE_H
