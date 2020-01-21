#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( application_suite )

UT_TEST_CASE( test1 )
{

  CKF_Triangulation topology( {3,3,3} );

  Visualizer vis;

  vis.add_topology(topology);

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( application_suite )
