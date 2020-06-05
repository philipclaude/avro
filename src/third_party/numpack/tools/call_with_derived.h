// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_CALL_WITH_DERIVED_H
#define SANS_CALL_WITH_DERIVED_H

#include <type_traits>
#include <typeinfo>
#include <utility> //std::forward

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
template < class iter, class end, class Func, class Base>
struct call_with_derived
{
  typedef typename boost::mpl::deref<iter>::type currentType;
  typedef typename boost::mpl::next< iter >::type next;

  template< class... Args >
  static bool call( Func&& f, const std::type_info& type, Base& base, Args&&... args )
  {

    // If the current type matches the derived type, then call f with the base class casted
    if (typeid(currentType) == type)
    {
      f(base.template cast<currentType>(), std::forward<Args>(args)...);
      return false;
    }

    // Check the next type to see if it's the derived type
    return call_with_derived< next, end, Func, Base >::call(std::forward<Func>(f), type, base, std::forward<Args>(args)...);
  }
};

// Recursive loop terminates when 'iter' is the same as 'end'
template < class end, class Func, class Base >
struct call_with_derived< end, end, Func, Base >
{
  template<class... Args>
  static bool call( Func&& f, const std::type_info& type, Base& base, Args&&... args )
  {
    // Did not find the derived type of base in the sequence
    return true;
  }
};

} // namespace detail


template < class Sequence, class Func, class Base, class... Args >
void call_with_derived( Func&& f, Base& base, Args&&... args )
{
  // Get the start and end iterators from the sequence
  typedef typename boost::mpl::begin< Sequence >::type begin;
  typedef typename boost::mpl::end< Sequence >::type   end;

  // Call the implementation and throw an exception if the derived class is not in Sequence
  if ( detail::call_with_derived< begin, end, Func, Base >::call(std::forward<Func>(f), base.derivedTypeID(), base, std::forward<Args>(args)...) )
    SANS_DEVELOPER_EXCEPTION("\nCould not find derived type:\n\n%s\n\nin the sequence:\n\n%s",
                             demangle(base.derivedTypeID().name()).c_str(), demangle(typeid(Sequence).name()).c_str() );
}

#define FRIEND_CALL_WITH_DERIVED \
template < class iter, class end, class Func, class Base> \
friend struct detail::call_with_derived;

}

#endif //SANS_CALL_WITH_DERIVED_H
