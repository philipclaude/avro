#include "unit_tester.hpp"

#include "graphics/raytracer.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_raytracer_suite )

UT_TEST_CASE( test1 )
{
  coord_t number = 3;
  coord_t dim = number;

  EGADS::Context context;
  std::vector<real_t> lengths(number,1.0);
  EGADS::Cube geometry(&context,lengths);

  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology( dims );
  topology.points().attach(geometry);

  int width = 500, height = width;

  RayTracer raytracer(width,height);

  raytracer.window().compute_view();
  raytracer.window().enable_controls(true);

  // initial draw, subsequent drawing will only be performed when a callback is invoked
  raytracer.draw();
  while (true) {

    // draw the scene
    raytracer.draw();

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(raytracer.window().window())) break;
    if (glfwGetKey(raytracer.window().window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( graphics_raytracer_suite )
