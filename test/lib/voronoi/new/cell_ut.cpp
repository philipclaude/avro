#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/new/cell.h"
#include "voronoi/new/diagram.h"

#define SPHERE 1

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_2d )
{
  GEO::PCK::initialize();

  static coord_t number = 2;
  static coord_t dim = number;
  index_t nb_points = 2e3;

  index_t N = 100;
  std::vector<index_t> dims(number,N);
  CKF_Triangulation ckf(dims);

  #if SPHERE
  Points sphere_points(3);
  real_t p[3];
  for (index_t k = 0; k < ckf.points().nb(); k++) {

    real_t theta = ckf.points()[k][0]*2.0*M_PI;
    real_t phi = ckf.points()[k][1]*M_PI;

    p[0] = cos(theta)*sin(phi);
    p[1] = sin(theta)*sin(phi);
    p[2] = cos(phi);
    sphere_points.create(p);
  }
  ckf.points().clear();
  ckf.points().set_dim(3);
  sphere_points.copy(ckf.points());

  // stitch together the 0 and 2pi boundaries
  for (index_t k = 0; k < N-1; k++) {

    // TODO

  }

  dim = 3;
  #endif

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

  Points sites(dim);

  #if 0
  ckf.points().copy( sites );
  nb_points = sites.nb();
  #elif SPHERE == 0
  for (index_t k = 0; k < nb_points; k++) {
    for (coord_t d = 0; d < dim; d++)
      x0[d] = random_within(0.0,1.0);
    sites.create(x0);
  }
  #else
  for (index_t k = 0; k < nb_points; k++) {
    real_t theta = random_within(0.0,2.0*M_PI);
    real_t phi = random_within(0.0,M_PI);

    p[0] = cos(theta)*sin(phi);
    p[1] = sin(theta)*sin(phi);
    p[2] = cos(phi);

    sites.create(p);
  }
  #endif

  voronoi::PowerDiagram diagram(topology,dim);

  #if SPHERE == 1
  diagram.set_ambient_dimension(3);
  #endif

  diagram.set_sites( sites );
  diagram.initialize();

  printf("computing cells\n");
  diagram.compute();
  printf("done computing cells\n");

  if (nb_points > 3e5 && number > 2) return;
  diagram.accumulate();
  diagram.create_field();

  printf("volume = %g\n",diagram.volume());
  printf("boundary area = %g\n",diagram.boundary_area());

  graphics::Viewer vis;
  vis.add(diagram);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_cell_test_suite )
