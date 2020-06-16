// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "common/types.h"

#include <numpack/types/SurrealD.h>
#include <numpack/types/SurrealS.h>

template<typename type>
type
func( type *x )
{
	return x[0]*x[1] +x[2]*sin(x[3] +x[4]*x[5]);
}

template<typename type>
type
func2(type* x )
{
	return x[0]*x[0]*x[1] +x[1]*x[1]*x[1];
}

using namespace avro;

UT_TEST_SUITE(surreal_test_suite)

UT_TEST_CASE(basic)
{
	std::vector< SurrealS<6> > x(6);
	real_t vals[6] = {2.,1.,3.4,2.3,0.1,2.10213};
	real_t ders[6] = {0,0,0,0,0,0};
	for (index_t k=0;k<6;k++)
	{
		ders[k] = 1.;
		x[k] = SurrealS<6>( vals[k] , ders , 6 );
		ders[k] = 0.;
	}
	SurrealS<6> f = func( x.data() );
	std::cout << f << std::endl;

	// simple test
	std::vector<SurrealS<2>> y(2);
	vals[0] = 1.;
	vals[1] = 2.;
	ders[0] = 1.;
	y[0] = SurrealS<2>( vals[0] , ders , 2);
	ders[0] = 0.;
	ders[1] = 1.;
	y[1] = SurrealS<2>( vals[1] , ders , 2);
	ders[1] = 0.;

	SurrealS<2> g = func2(y.data());
	std::cout << g << std::endl;

	printf("df/dx = %g, df/dy = %g\n",g.deriv(0),g.deriv(1));
}
UT_TEST_CASE_END(basic)

UT_TEST_SUITE_END(surreal_test_suite)
