#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/cell.h"
#include "voronoi/diagram.h"

#define SPHERE 0

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_2d )
{
  GEO::PCK::initialize();

  static coord_t number = 2;
  static coord_t dim = number+1;
  index_t nb_points = 1e4;

  index_t N = 50;
  std::vector<index_t> dims(number,N);
  CKF_Triangulation ckf(dims);

  real_t point[dim+1];
  for (coord_t d = 0; d < dim+1; d++)
    point[d] = 0.0;

  #if SPHERE
  // stitch together the 0 and 2pi boundaries
  std::map<index_t,index_t> point_map;
  for (index_t k = 0; k < N; k++) {
    index_t p = k;
    index_t q = (N-1)*N + k;

    point_map.insert( {q,p} );

    p = N*k;
    point_map.insert( {p,0} );

    p += N-1;
    point_map.insert( {p,N-1} );
  }

  for (index_t k = 0; k < ckf.nb(); k++) {
    for (index_t j = 0; j < ckf.nv(k); j++) {
      index_t q = ckf(k,j);
      if (point_map.find(q) == point_map.end()) continue;
      ckf(k,j) = point_map.at(q);
    }
  }

  std::vector<real_t> volumes;
  ckf.get_volumes(volumes);
  std::vector<index_t> removals;
  for (index_t k = 0; k < ckf.nb(); k++) {
    if (fabs(volumes[k]) < 1e-8) {
      removals.push_back(k);
    }
  }

  std::sort( removals.begin() , removals.end () );
  for (index_t k = 0; k < removals.size(); k++) {
    ckf.remove( removals[k] - k );
  }
  ckf.remove_unused();

  Points sphere_points(dim+1);
  for (index_t k = 0; k < ckf.points().nb(); k++) {

    real_t theta = ckf.points()[k][0]*2.0*M_PI;
    real_t phi = ckf.points()[k][1]*M_PI;

    point[0] = cos(theta)*sin(phi);
    point[1] = sin(theta)*sin(phi);
    point[2] = cos(phi);
    sphere_points.create(point);
  }
  ckf.points().clear();
  ckf.points().set_dim(dim+1);
  sphere_points.copy(ckf.points());

  dim += 1;
  #endif

  // the domain which will be used to clip the voronoi diagram
  Topology<Simplex>& topology = ckf;

  Points sites(dim);
  #if 0 // test exactness
  ckf.points().copy( sites );
  nb_points = sites.nb();
  #elif SPHERE == 0
  for (index_t k = 0; k < nb_points; k++) {
    for (coord_t d = 0; d < number; d++)
      point[d] = random_within(0.0,1.0);
    sites.create(point);
  }
  #else
  for (index_t k = 0; k < nb_points; k++) {

    // compute a random point in the parameter space of the sphere
    real_t theta = random_within(0.0,2.0*M_PI);
    real_t phi = random_within(0.0,M_PI);

    // map the point to the sphere
    point[0] = cos(theta)*sin(phi);
    point[1] = sin(theta)*sin(phi);
    point[2] = cos(phi);

    // create the delaunay site
    sites.create(point);
  }

  #endif

  voronoi::PowerDiagram diagram(topology,dim);

  #if SPHERE == 1
  diagram.set_ambient_dimension(3);
  #endif

  diagram.set_sites( sites );
  diagram.initialize();

  diagram.optimize_points(200);

  std::vector<real_t> mass( sites.nb() , diagram.volume() / real_t(sites.nb()) );

  //mass = diagram.cell_volume();
  diagram.optimize_weights_kmt( 20 , mass );

  diagram.accumulate();
  diagram.create_field();

  real_t volume = diagram.volume();
  real_t area = diagram.area();
  printf("volume = %g\n",volume);
  printf("area   = %g\n",area);

  #if SPHERE == 0
  UT_ASSERT_NEAR( volume , 1.0 , 1e-12 );
  UT_ASSERT_NEAR( area , 2.0*number , 1e-12 );
  #else
  UT_ASSERT_NEAR( volume , 4.*M_PI , 1e-12 );
  UT_ASSERT_NEAR( area , 0.0 , 1e-12 );
  #endif

  graphics::Viewer vis;
  vis.add(diagram);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_cell_test_suite )
