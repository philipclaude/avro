#include "unit_tester.hpp"

#include "graphics/raytracer.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/factory.h"

#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_raytracer_suite )

vec3
get_center(const Points& points) {
  index_t np = points.nb();
  vec3 c = {0.,0.,0.};
  for (index_t i = 0; i < np; i++) {
    for (coord_t d = 0; d < 3; d++)
      c[d] += points(i,d) / np;
  }
  //printf("center = \n");
  //c.print();
  return c;
}

UT_TEST_CASE( test1 )
{

  // get the input mesh
  std::string meshname( "library/meshes/giraffe.avro" );
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh = library::get_mesh(meshname,ptopology);
  Topology<Simplex>& topology = *static_cast<Topology<Simplex>*>(ptopology.get());

  int width = 600, height = 0.6*width;

  RayTracer raytracer(width,height);

  raytracer.window().camera().set_lookat( get_center( topology.points() ) );

  // first sphere: red ball
  Material red;
  red.type = 2; // refractive
  red.eta  = 1.5;
  red.ka = {0.4,0.,0.};
  red.kd = {0.4,0.,0.};
  red.ks = {0.4,0.,0.};
  Sphere sphere_b( {4.,1.,2.} , 1.0 , red );
  raytracer.scene().add( sphere_b );

  // second sphere: ground
  Material ground;
  ground.ka = {0.25,0.25,0.25};
  ground.kd = {0.25,0.25,0.25};
  ground.ks = {0.25,0.25,0.25};
  ground.shine = 0.0;
  Sphere sphere_g( {0.,-1000.,0.} , 1000. , ground );
  raytracer.scene().add( sphere_g );

  // third sphere: mint ball
  Material mater2;
  mater2.type = 0; // diffuse
  mater2.ka = {0.07,0.98,0.53};
  mater2.kd = mater2.ka;
  mater2.ks = mater2.kd;
  Sphere sphere2( {-1.,0.5,2.} , 0.5 , mater2 );
  raytracer.scene().add( sphere2 );

  // fourth sphere: salmon ball
  Material mater3;
  mater3.type = 0; // diffuse
  mater3.ka = {0.98,0.5,0.44};
  mater3.kd = mater3.ka;
  mater3.ks = mater3.kd;
  Sphere sphere3( {1.,0.5,-0.5} , 0.5 , mater3 );
  raytracer.scene().add( sphere3 );

  // mesh
  Material mater4;
  mater4.type = 0; // reflective
  mater4.ka = {0.98,0.5,0.44};
  mater4.kd = mater4.ka;
  mater4.ks = mater4.kd;
  raytracer.scene().add( topology , mater4 , true ); // true means use BVH

  // point light source (white)
  Light light;
  light.La = {1.,1.,1.};
  light.Ld = {1.,1.,1.};
  light.Ls = {1.,1.,1.};
  light.position = {-10.,5.,0.};
  raytracer.lights().push_back( light );

  raytracer.window().enable_controls(true);

  // initial draw, subsequent drawing will only be performed when a callback is invoked
  raytracer.render();
  while (true) {

    // draw the scene
    //raytracer.render();

    // wait for user input
    glfwWaitEvents();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(raytracer.window().window())) break;
    if (glfwGetKey(raytracer.window().window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( graphics_raytracer_suite )
