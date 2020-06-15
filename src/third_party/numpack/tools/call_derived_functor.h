// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_CALL_DERIVED_FUNCTOR_H
#define SANS_CALL_DERIVED_FUNCTOR_H

#include <type_traits>
#include <typeinfo>
#include <utility>

#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/deref.hpp>

#include "SANSException.h"
#include "demangle.h"

namespace numpack 
{

namespace detail
{

// Looping class. Check the current derived type stored in 'iter' to see if the 'base' is a base class of it.
// Otherwise continue to the next derived class in the sequence.
template < class iter, class end, class Base>
struct call_derived_functor
{
  static bool getIndex( Base& base, int &index )
  {
    typedef typename boost::mpl::deref<iter>::type currentType;

    // If the current type matches the derived type, then the index is already set and simply return a success
    if (typeid(currentType) == base.derivedTypeID())
    {
      return false;
    }

    // Increase the index
    index++;

    // Check the next type to see if it's the derived type
    return call_derived_functor< typename boost::mpl::next< iter >::type, end, Base >::getIndex(base, index);
  }

  template<class ...Args>
  static bool call( Base& base, const int &index, int& count, Args &&... args )
  {
    typedef typename boost::mpl::deref<iter>::type currentType;

    // If the index matches the current count, then we have the right derived type
    if (index == count)
    {
      base.template cast<currentType>().operator()(std::forward<Args>(args)...);
      return false;
    }

    // Increase the count
    count++;

    // Check the next type to see if it's the right count
    return call_derived_functor< typename boost::mpl::next< iter >::type, end, Base >::template call(base, index, count, std::forward<Args>(args)...);
  }
};

// Recursive loop terminates when 'iter' is the same as 'end'
template < class end, class Base >
struct call_derived_functor< end, end, Base >
{
  static bool getIndex( Base& base, int &index )
  {
    // Did not find the derived type of base in the sequence
    return true;
  }

  template<class ...Args>
  static bool call( Base& base, const int &index, int& count, Args &&... args )
  {
    // Did not find the derived type of base in the sequence
    return true;
  }
};

} // namespace detail


template < class Sequence, class Base >
struct call_derived_functor
{
  // Get the start and end iterators from the sequence
  typedef typename boost::mpl::begin< Sequence >::type begin;
  typedef typename boost::mpl::end< Sequence >::type   end;

  // cppcheck-suppress noExplicitConstructor
  call_derived_functor( Base& base ) : index_(0), base_(base)
  {
    // Find the index for the derived type and throw an exception if the derived class is not in Sequence
    if ( detail::call_derived_functor< begin, end, Base >::getIndex(base_, index_) )
      SANS_DEVELOPER_EXCEPTION("\nCould not find derived type:\n\n%s\n\nin the sequence:\n\n%s",
                               demangle(base.derivedTypeID().name()).c_str(), demangle(typeid(Sequence).name()).c_str() );
  }

  template<class ...Args>
  void operator()( Args &&... args ) const
  {
    int count = 0;
    // Call the implementation and throw an exception if the derived class is not in Sequence
    if ( detail::call_derived_functor< begin, end, Base >::template call( base_, index_, count, std::forward<Args>(args)...) )
      SANS_DEVELOPER_EXCEPTION("\nCould not find derived type:\n\n%s\n\nin the sequence:\n\n%s",
                               demangle(base_.derivedTypeID().name()).c_str(), demangle(typeid(Sequence).name()).c_str() );

  }

protected:
  int index_;
  Base& base_;
};


}

#endif //SANS_CALL_DERIVED_FUNCTOR_H
