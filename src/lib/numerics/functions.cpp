//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "numerics/functions.h"

#include <cmath> // nextafter

namespace avro
{

namespace numerics
{

index_t
factorial( const index_t i ) {
  if (i == 0 || i == 1) return 1;
  return i*factorial(i-1);
}

index_t
binomial( const index_t n , const index_t k ) {
	index_t num = 1,den = 1;
	for (index_t i = 1; i <= k; i++) {
		num *= (n +1 -i);
		den *= i;
	}
	return num/den;;
}

index_t
nchoosek( index_t n , index_t k ) {
  return factorial(n)/( factorial(k)*factorial(n-k) );
}

real_t
eps( const real_t& x ) {
  // this is like matlab's eps function
  return ::nextafter( x , x +1.0f ) -x;
}

} // numerics

} // avro
