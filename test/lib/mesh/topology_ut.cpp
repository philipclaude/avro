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

#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/tesseract.h"

#include "shape/shape.h"
#include "shape/quadrature.h"

#include "mesh/topology.h"
#include "mesh/points.h"

using namespace avro;

UT_TEST_SUITE( mesh_topology_suite )

UT_TEST_CASE( simplex_tests )
{
  Points vertices(3);
  coord_t number = 3;

  Topology<Simplex> topology(vertices,number);

  std::shared_ptr< Topology<Simplex> > leaf = std::make_shared<Topology<Simplex>>(vertices,number);
  topology.add_child(leaf);

  std::shared_ptr< Topology<Simplex> > node = std::make_shared<Topology<Simplex>>(vertices,0);
  leaf->add_child(node);

  Topology<Simplex>& c = topology.topology(0);
  UNUSED(c);

  Topology<Simplex> topology_copy( vertices , number );
  topology_copy.Tree<Topology<Simplex>>::copy(topology);
  
  topology_copy.Tree<Topology<Simplex>>::print();

  topology_copy.Tree<Topology<Simplex>>::print();


  UT_ASSERT_EQUALS( topology_copy.nb_children() , 1 );
  UT_ASSERT_EQUALS( topology_copy.child(0).nb_children() , 1 );
}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_CASE( hierarchy_from_geometry )
{
  EGADS::Context context;
  EGADS::Cube box( &context, {1,1,1} );

  Points points(3);
  Topology<Simplex> topology(points,2);
  topology.Tree<Topology<Simplex>>::copy( *box.child(0) );

  UT_ASSERT_EQUALS( box.child(0)->nb_children() , 6 );
  UT_ASSERT_EQUALS( topology.nb_children() , 6 );

  box.child(0)->Tree<Entity>::print();
  topology.Tree<Topology<Simplex>>::print();

  for (index_t k=0;k<box.child(0)->nb_children();k++)
  {
    UT_ASSERT_EQUALS( box.child(0)->child(k).nb_children() , 1 );
    UT_ASSERT_EQUALS( box.child(0)->child(k).child(0).nb_children() , 4 );

    UT_ASSERT_EQUALS( topology.child(k).nb_children() , 1 );
    UT_ASSERT_EQUALS( topology.child(k).child(0).nb_children() , 4 );
  }

}
UT_TEST_CASE_END( hierarchy_from_geometry )

UT_TEST_CASE( simplex_close )
{
  for (coord_t dim=2;dim<=4;dim++)
  {
    for (index_t N=4;N<=4;N+=2)
    {
      std::vector<index_t> dims(dim,N);
      CKF_Triangulation topology(dims);

      topology.close();

      Topology<Simplex> boundary( topology.points() , topology.number()-1 );
      topology.get_boundary(boundary);
      UT_ASSERT_EQUALS( boundary.nb() , 0 );
    }
  }

}
UT_TEST_CASE_END( simplex_close )

UT_TEST_SUITE_END( mesh_topology_suite )
