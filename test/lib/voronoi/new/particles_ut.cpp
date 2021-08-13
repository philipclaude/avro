#include "unit_tester.hpp"

#include "graphics/application.h"

#include "mesh/points.h"
#include "voronoi/new/cell.h"
#include "voronoi/new/particles.h"

using namespace avro;

UT_TEST_SUITE( voronoi_particles_test_suite )

UT_TEST_CASE( test_2d )
{

  index_t nb_particles = 1e2;

  voronoi::ParticleSimulator particles( "CKF-4-4" , nb_particles );

  particles.sample("random","lloyd");
  particles.set_density();

  particles.simulate(10);

  graphics::Viewer vis;
  vis.add(particles.diagram());

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_particles_test_suite )
