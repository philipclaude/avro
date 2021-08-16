#include "unit_tester.hpp"

#include "../bin/programs.h"

#include "graphics/application.h"

#include "mesh/points.h"
#include "voronoi/cell.h"
#include "voronoi/particles.h"

using namespace avro;

UT_TEST_SUITE( voronoi_particles_test_suite )

UT_TEST_CASE( test_2d )
{
  index_t nb_particles = 1e3;

  voronoi::ParticleSimulator particles( "CKF-2-2" , nb_particles );

  particles.sample("random","lloyd");
  particles.set_density();

  particles.save_every(20);
  particles.simulate(1200);
  //particles.save_frames("particles.json");

  graphics::Viewer vis;
  vis.add(particles.diagram());
  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_particles_test_suite )
