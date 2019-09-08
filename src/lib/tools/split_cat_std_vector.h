// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SPLIT_CAT_STD_VECTOR_H
#define SPLIT_CAT_STD_VECTOR_H

#include <vector>

namespace SANS
{

template< class T >
const std::vector<T>
cat( const std::vector<std::vector<T>>& v )
{
  std::vector<T> vec ={};
  for (auto it = v.begin(); it != v.end(); ++it)
    vec.insert( vec.end(), it->begin(), it->end() );
  return vec;
}

template< class T >
const std::vector<std::vector<T>>
split( const std::vector<T>& v )
{
  std::vector<std::vector<T>> vec;
  for (auto it = v.begin(); it != v.end(); ++it)
    vec.push_back(std::vector<T>(1,*it));
  return vec;
}

}
#endif // SPLIT_CAT_STD_VECTOR_H
