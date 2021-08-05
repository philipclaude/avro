#include "unit_tester.hpp"

#include "graphics/bsp.h"
#include "graphics/managers.h"
#include "graphics/vao.h"
#include "graphics/window.h"

#include "graphics/colormap.h"
#include "graphics/shader.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_bsp_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 3;

  EGADS::Context context;
  std::vector<real_t> lengths(number,1.0);
  EGADS::Cube geometry(&context,lengths);

  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology( dims );
  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  topology.points().attach(geometry);

  if (AVRO_FULL_UNIT_TEST) return;

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

  BSPTree tree;
  tree.build(triangles);

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( graphics_bsp_suite )
