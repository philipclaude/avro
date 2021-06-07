//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "common/types.h"

#include "mesh/points.h"

#include "numerics/geometry.h"

#if 0
extern "C"
{
#include <qhull/libqhull/libqhull.h>
#include <qhull/libqhull/geom.h>
}
#endif

using namespace avro;


UT_TEST_SUITE( qhull_test_suite )

UT_TEST_CASE( test_nd )
{
  #if 0
  coord_t dim = 2;

  for (index_t np = 3; np <= 10; np++)
  {
    real_t  dt = 2.*M_PI/(np);

    Points points(dim);
    for (index_t k=0;k<np;k++)
    {
      real_t x[2] = { cos(k*dt) , sin(k*dt) };
      points.create(x);
    }

    real_t s = numerics::distance( points[0] , points[1] , dim );
    printf("s = %g\n",s);

    points.print();

    std::vector<coordT> coordinates(np*dim);

    index_t i = 0;
    for (index_t k=0;k<points.nb();k++)
    for (coord_t d=0;d<dim;d++)
      coordinates[i++] = points[k][d];


    char flags[25];
    sprintf (flags, "qhull s FA");

    qh_new_qhull(dim, np , coordinates.data() , 0, flags, NULL, NULL);
    qh_getarea(qh facet_list);
    real_t area  = (qh totvol);
    real_t perim = (qh totarea);
    qh_freeqhull(!qh_ALL);

    printf("area = %g, expect: %g\n",area,0.25*s*s*np/tan(M_PI/np));
    printf("perim = %g, expect: %g\n",perim,np*s);

    UT_ASSERT_NEAR( area , 0.25*s*s*np/tan(M_PI/np) , 1e-12 );
    UT_ASSERT_NEAR( perim , np*s , 1e-12 );
  }

  #endif

}
UT_TEST_CASE_END( test_nd )

UT_TEST_SUITE_END( qhull_test_suite )
