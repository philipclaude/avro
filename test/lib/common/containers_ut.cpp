//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "common/tools.h"
#include "avro_types.h"

#include <list>
#include <string>
#include <vector>
#include <time.h>

using namespace avro;

UT_TEST_SUITE(containers_test_suite)


UT_TEST_CASE(test1)
{
  clock_t t0, t1;
  index_t nb_elem = 4e5;
  std::vector<index_t> data(nb_elem,0);

  t0 = clock();
  data.resize( nb_elem , 0 );
  std::list<index_t> ldata(data.begin(),data.end());
  for (index_t k = 0; k < nb_elem-1; k++)
    ldata.erase( ldata.begin() );
  t1 = clock();
  printf("list erase time: %3.4g\n",real_t(t1-t0)/real_t(CLOCKS_PER_SEC));


  data.resize( nb_elem , 0 );
  t0 = clock();
  for (index_t k = 0; k < nb_elem-1; k++) {
    data.erase( data.begin() );
  }
  t1 = clock();
  printf("vector erase time: %3.4g\n",real_t(t1-t0)/real_t(CLOCKS_PER_SEC));

  data.resize( nb_elem , 0 );


}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(containers_test_suite)
