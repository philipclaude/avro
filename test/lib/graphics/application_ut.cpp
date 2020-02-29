#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/obj.h"
#include "library/samples.h"

#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( application_suite )

UT_TEST_CASE( test1 )
{

  #if 1
  library::objFile topology( "/Users/pcaplan/Google Drive/library/models/obj/spot.obj" );
  #elif 0
  CKF_Triangulation topology( {4,4,4} );
  #endif

  Delaunay delaunay(topology.points().dim());
  topology.points().copy(delaunay);

  Delaunay z(topology.points());
  for (index_t k=0;k<z.nb();k++)
  for (index_t d=0;d<z.dim();d++)
    z[k][d] = random_within( 0.0 , 1.0 );

  delaunay::RestrictedVoronoiDiagram rvd(topology,z);
  rvd.compute(true);

  Visualizer vis;

  vis.add_topology(topology);
  //vis.add_topology(rvd);

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( application_suite )
