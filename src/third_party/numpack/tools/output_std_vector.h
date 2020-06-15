// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef OUTPUT_STD_VECTOR_H
#define OUTPUT_STD_VECTOR_H

#include "tools/SANSnumerics.h"     // Real

#include <iostream>
#include <vector>

namespace numpack 
{

template< class T >
std::ostream&
operator<<( std::ostream& os, const std::vector<T>& v )
{
  const std::size_t size = v.size();
  for (std::size_t ii = 0; ii < size; ii++)
  {
    os << v[ii];
    if (ii < size-1) os << ", ";
  }
  return os;
}

template< class T >
std::ostream&
operator<<( std::ostream& os, const std::vector<std::vector<T>>& v )
{
  for (std::size_t ii = 0; ii < v.size(); ii++)
  {
    const std::size_t size = v[ii].size();
    for (std::size_t jj = 0; jj < size; jj++)
    {
      os << v[ii][jj];
      if (jj < size-1) os << ", ";
    }
    if (ii < v.size()-1) os << ", ";
  }
  return os;
}

}
#endif // OUTPUT_STD_VECTOR_H
