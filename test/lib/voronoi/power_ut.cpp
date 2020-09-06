#include "unit_tester.hpp"

#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/delaunay.h"
#include "voronoi/power.h"
#include "voronoi/voronoi_cell.h"

using namespace avro;

UT_TEST_SUITE( laguerre_test_suite )

UT_TEST_CASE( test_simplex )
{
  GEO::PCK::initialize();

  coord_t number = 2;
  coord_t dim = number;

  std::vector<index_t> sizes(number,2);
  CKF_Triangulation domain(sizes);

  // create random delaunay vertices
  Delaunay delaunay( dim );
  #if 1
  index_t nb_points = 1e2;
  std::vector<real_t> x(dim,0.);
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t d=0;d<dim;d++)
      x[d] = random_within(0.,1.);
    delaunay.create(x.data());
  }
  #else
  std::vector<index_t> dims(number,10);
  CKF_Triangulation ckf(dims);
  ckf.points().copy(delaunay);
  #endif

  std::vector<real_t> weights( delaunay.nb() , 0.0 );
  delaunay::PowerDiagram diagram( delaunay , domain , weights );
  diagram.set_exact(true);
  diagram.compute();

  diagram.optimize_cvt();

  printf("optimizing transport map:\n");
  diagram.optimize_otm();

  graphics::Visualizer vis;
  vis.add_topology(diagram);

  if (number<4)
    vis.run();
}
UT_TEST_CASE_END( test_simplex )

UT_TEST_SUITE_END( laguerre_test_suite )
