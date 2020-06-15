// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_SAFE_AT_H
#define SANS_SAFE_AT_H

#include "SANSException.h"

#include <sstream>
#include <type_traits>
#include <stdexcept> //std::out_of_range
#include <boost/preprocessor/stringize.hpp>

namespace numpack 
{

// Always use this macro so file and line numbers get included in the error message
#define safe_at( container, key) safe_at_(__FILE__,__LINE__, BOOST_PP_STRINGIZE(container), BOOST_PP_STRINGIZE(key), container, key )

// std::map
template< template<class,class,class,class> class container_type, class Key, class T, class Compare, class Allocator>
inline typename std::conditional<std::is_const<container_type<Key, T, Compare, Allocator>>::value, const T&, T&>::type
safe_at_(const char* file, const int line, const char* container_str, const char* key_str,
         container_type<Key, T, Compare, Allocator>& container, const Key& key)
{
  try
  {
    return container.at(key);
  }
  catch ( const std::out_of_range& e )
  {
    std::stringstream msg;
    msg << std::endl;
    msg << container_str << ".at(" << key_str << ") does not have key == " << key << std::endl;
    msg << std::endl;
    msg << file << ": " << line;
    SANS_DEVELOPER_EXCEPTION( msg.str() );
  }
}

// std::unordered_map
template< template<class,class,class,class,class> class container_type, class Key, class T, class Hash, class Compare, class Allocator>
inline typename std::conditional<std::is_const<container_type<Key, T, Hash, Compare, Allocator>>::value, const T&, T&>::type
safe_at_(const char* file, const int line, const char* container_str, const char* key_str,
         container_type<Key, T, Hash, Compare, Allocator>& container, const Key& key)
{
  try
  {
    return container.at(key);
  }
  catch ( const std::out_of_range& e )
  {
    std::stringstream msg;
    msg << std::endl;
    msg << container_str << ".at(" << key_str << ") does not have key == " << key << std::endl;
    msg << std::endl;
    msg << file << ": " << line;
    SANS_DEVELOPER_EXCEPTION( msg.str() );
  }
}

}

#endif //SANS_SAFE_AT_H
