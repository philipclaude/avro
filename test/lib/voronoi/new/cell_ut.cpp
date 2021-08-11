#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/field.hpp"

#include "voronoi/new/cell.h"
#include "voronoi/new/diagram.h"


using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_2d )
{
  GEO::PCK::initialize();

  static coord_t number = 2;
  static coord_t dim = number;
  index_t nb_points = 1e2;

  std::vector<index_t> dims(number,10);
  CKF_Triangulation ckf(dims);

  Points points(dim);
  real_t x0[] = {0,0,0}; points.create(x0);
  real_t x1[] = {1,0,0}; points.create(x1);
  real_t x2[] = {0,1,0}; points.create(x2);
  if (number > 2) {
    real_t x3[] = {0,0,1}; points.create(x3);
  }

  #if 0
  Topology<Simplex> topology(points,number);
  index_t t[] = {0,1,2,3};
  topology.add(t,number+1);
  #else
  Topology<Simplex>& topology = ckf;
  #endif
  topology.build_structures();

  Points sites(dim);

  #if 1
  ckf.points().copy( sites );
  nb_points = sites.nb();
  #else
  for (index_t k = 0; k < nb_points; k++) {
    for (coord_t d = 0; d < dim; d++)
      x0[d] = random_within(0.0,1.0);
    sites.create(x0);
  }
  #endif

  voronoi::PowerDiagram diagram(topology,dim);
  diagram.set_sites( sites );
  diagram.initialize( sites.nb() );


  printf("computing cells\n");
  diagram.compute();
  printf("done computing cells\n");
  diagram.accumulate();
  diagram.create_field();

  printf("volume = %g\n",diagram.volume());

  graphics::Viewer vis;
  vis.add(diagram);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_cell_test_suite )
