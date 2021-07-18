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

#include "geometry/entity.h"
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/csm.h"

using namespace avro;

UT_TEST_SUITE( csm_test_suite )

UT_TEST_CASE(test1)
{

  OpenCSM_Model model(BASE_TEST_DIR+"/geometry/bottle.csm");

  TessellationParameters params;
  params.set("min size" , 0.1 );
  params.set("min angle", 20.0);

  ModelTessellation tess(model,params);

  for (index_t k=0;k<tess.points().nb();k++)
  {
    if (tess.points().entity(k)->number()<=1)
      UT_ASSERT( tess.points().u(k,1) > 1e10 );
    else
      UT_ASSERT( tess.points().u(k,1) < 1e10 );
  }

#if 1
  graphics::Viewer vis;
  vis.add( tess.topology(0) );

#else
  graphics::WebVisualizer vis;

  Topology<Simplex>& topology = static_cast<Topology<Simplex>&>( tess.topology(0) );

  printf("topology nb = %lu\n",topology.nb());
  topology.Tree<Topology<Simplex>>::print();
  tess.topology(0).Table<index_t>::print();

  for (index_t k=0;k<15;k++)//topology.nb_children();k++)
    vis.add( topology.child(k) );

#endif

  vis.run();

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END( csm_test_suite )
