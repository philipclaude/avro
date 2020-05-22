#include "unit_tester.hpp"

#include "graphics/primitive.h"

#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;

UT_TEST_SUITE( graphics_primitive_suite )

UT_TEST_CASE( tree_test )
{
  Points vertices(3);
  coord_t number = 3;

  Topology<Simplex> topology(vertices,number);

  std::shared_ptr< Topology<Simplex> > leaf = std::make_shared<Topology<Simplex>>(vertices,number);
  topology.add_child(leaf);

  std::shared_ptr< Topology<Simplex> > node = std::make_shared<Topology<Simplex>>(vertices,0);
  leaf->add_child(node);

  graphics::Primitive primitive( topology , nullptr );
  primitive.copy( topology );

  UT_ASSERT_EQUALS( primitive.nb_children() , 1 );
  UT_ASSERT_EQUALS( primitive.child(0).nb_children() , 1 );
}
UT_TEST_CASE_END( tree_test )

UT_TEST_SUITE_END( graphics_primitive_suite )
