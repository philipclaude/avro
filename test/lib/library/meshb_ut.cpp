#include "unit_tester.hpp"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/meshb.h"

using namespace avro;

UT_TEST_SUITE( meshb_test_suite )

UT_TEST_CASE( test1 )
{
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");

  library::meshb mesh( BASE_TEST_DIR+"/meshes/cube-cylinder.mesh" , &model );
  mesh.points().compute_params();

  mesh.write( mesh , "tmp/cc.mesh" , true );
  library::meshb mesh_in( "tmp/cc.mesh" );
  UT_ASSERT_EQUALS( mesh_in.nb_topologies() , 8 ); // 7 faces + 1 volume

  graphics::Visualizer vis;

  std::shared_ptr<graphics::Widget> toolbar = std::make_shared<graphics::Toolbar>(vis.main_window(),vis);
  vis.main_window().interface().add_widget( toolbar );

  for (index_t k=0;k<mesh_in.nb_topologies();k++)
    vis.add_topology(mesh_in.topology(k));

  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( meshb_test_suite )
