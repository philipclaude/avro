//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_COMMON_SET_H_
#define AVRO_COMMON_SET_H_

#include "common/tools.h"
#include "common/types.h"

#include <vector>
#include <algorithm>

namespace avro
{

namespace Set
{

template <typename T>
int
contains( const T& x , const std::vector<T>& X )
{
	const index_t n = X.size();
	for (index_t k=0;k<n;k++)
		if (X[k]==x) return int(k);
	return -1;
}

template<typename T>
void
intersection( const std::vector<T>& s0 , const std::vector<T>& s1 , std::vector<T>& p )
{
	p.clear();
	const index_t n0 = s0.size();

	p.clear();
	for (index_t i=0;i<n0;i++)
	{
		if ( contains( s0[i] , s1 ) > -1 )
			p.push_back(s0[i]);
	}
	uniquify(p);
}

} // Set

} // avro

#endif
