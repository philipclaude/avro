#include "unit_tester.hpp"

#include "graphics/bsp.h"
#include "graphics/managers.h"
#include "graphics/postscript.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/factory.h"

#include "mesh/field.hpp"
#include "mesh/mesh.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include <time.h>

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_bsp_suite )

UT_TEST_CASE( test1 )
{
  typedef Simplex type;
  std::string mesh_name = "CKF-4-4-4";//AVRO_SOURCE_DIR + "/build/release/ccp2_9.mesh";

  // get the original input mesh
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(mesh_name,ptopology);
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());

  //topology.points()[0][0] += 0.1;

  topology.orient();

  int width = 400, height = width;
  Window window(width,height);
  window.init();

  Plot plot(topology);
  plot.build();
  window.add_plot(&plot);

  window.compute_view();
  window.enable_controls(true);

  const mat4& view_matrix = window.camera().view_matrix();
  const mat4& projection_matrix = window.camera().projection_matrix();
  const mat4& screen_matrix = window.screen_matrix();

  screen_matrix.print();

  BSPTriangles triangles;

  triangles.build( plot , view_matrix , projection_matrix , screen_matrix );

  clock_t t0 = clock();
  BSPTree tree;
  tree.build(triangles);
  clock_t t1 = clock();
  printf("BSP time = %g sec.\n",real_t(t1-t0)/real_t(CLOCKS_PER_SEC));

  std::vector<BSPTriangle*> T;

  const vec3& eye = window.camera().eye();
  tree.get_triangles( eye , T );

  printf("retrieved %lu triangles\n",T.size());
  PostScriptWriter writer(AVRO_SOURCE_DIR + "/test/tmp/test.eps");
  writer.begin( window.width() , window.height() );
  writer.write( T , view_matrix , projection_matrix );
  writer.end();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( graphics_bsp_suite )
