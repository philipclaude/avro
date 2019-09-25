// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_ARRAY_H_
#define URSA_COMMON_ARRAY_H_

#include "common/error.h"
#include "common/types.h"
#include "common/stringify.h"

#include <vector>
#include <algorithm>

namespace ursa
{

template<typename type>
index_t
indexof( const type* x , const std::vector<type*>& X )
{
	for (index_t k=0;k<X.size();k++)
	{
		if (X[k]==x) return k;
	}
	ursa_assert_not_reached;
	return 0; // to avoid compiler warnings
}

template <typename T>
void
uniquify( std::vector<T> &x )
{
  // uniquify the list of indices
  std::sort( x.begin() , x.end() );
  typename std::vector<T>::iterator it = std::unique( x.begin() , x.end() );
  x.erase(it,x.end());
}

template <typename T>
index_t
uniquify_noerase( std::vector<T> &x , const index_t upto )
{
  // uniquify the list of indices
  std::sort( x.begin() , x.end() );
  typename std::vector<T>::iterator it = std::unique( x.begin() , x.begin()+upto );
  return std::distance(x.begin(),it);
}

template<typename T>
void
reverse( std::vector<T>& x )
{
	std::reverse( x.begin() , x.end() );
}

template<typename T>
std::vector<T>
linspace( T lo , T hi , int n=-1 )
{
  if (n==-1) n = hi -lo +1;
  T dt = T(hi -lo)/T(n -1);
  std::vector<T> y(n);
  for (index_t k=0;k<index_t(n);k++)
    y[k] = lo +k*dt;
  return y;
}

inline
std::vector<index_t>
linspace( index_t n )
{
	std::vector<index_t> y(n,0);
	for (index_t k=1;k<n;k++)
		y[k] = y[k-1] +1;
	return y;
}

template<typename T>
std::string
unique_label( T& x0 , T& x1 )
{
  ursa_assert(x0!=x1);
  if (x0>x1) std::swap(x0,x1);
  std::string s = stringify(x0)+"|"+stringify(x1);
  return s;
}

template<typename T>
std::string
unique_label( std::vector<T>& x )
{
  ursa_assert( x.size()>0 );
  std::sort( x.begin() , x.end() );
  std::string s = stringify(x[0]);
  for (index_t k=1;k<x.size();k++)
    s += "|"+stringify(x[k]);
  return s;
}

template<typename T>
std::string
unique_label_skip( T* x , index_t n , index_t kskip )
{
  std::vector<index_t> f(x,x+n);
  f.erase( f.begin() +kskip );
  return unique_label(f);
}

} // ursa

#endif
