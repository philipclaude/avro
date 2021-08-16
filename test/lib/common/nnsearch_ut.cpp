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

#include "common/tools.h"

#include <geogram/nn_search.h>

using namespace avro;

UT_TEST_SUITE(nnsearch_test_suite)

UT_TEST_CASE(test1)
{
  const coord_t dim = 4;
  GEO::NearestNeighborSearch* nns = GEO::NearestNeighborSearch::create(dim,"BNN");

  //nns->set_exact(false);

  index_t nb_points = 1e5;
  std::vector<real_t> x( nb_points*dim );
  for (index_t k=0;k<nb_points*dim;k++)
    x[k] = random_within(0.,1.);

  clock_t t0 = clock();
  nns->set_points( nb_points , x.data() );
  printf("time to set points = %g\n",real_t(clock()-t0)/CLOCKS_PER_SEC);

  index_t nb_neighbors = 50;
  std::vector<index_t> neighbors(nb_neighbors);
  std::vector<double> neighbors_sq_dist(nb_neighbors);

  //#pragma omp parallel for
  t0 = clock();
  for (index_t k=0;k<nb_points;k++)
  {
    nns->get_nearest_neighbors( nb_neighbors , k , neighbors.data() , neighbors_sq_dist.data() );
    nns->get_nearest_neighbors( nb_neighbors , &x[k*dim] , neighbors.data() , neighbors_sq_dist.data() , GEO::NearestNeighborSearch::KeepInitialValues() );
    UT_ASSERT_EQUALS( neighbors[0] , k );
    //print_inline(neighbors);
  }
  printf("time to compute neighbours = %g\n",real_t(clock()-t0)/CLOCKS_PER_SEC);

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END( nnsearch_test_suite )
