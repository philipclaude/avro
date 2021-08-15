#include "programs.h"

#include "common/process.h"
#include "common/tools.h"

#include "graphics/gl.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include "library/ckf.h"

#include <json/json.hpp>

#include <fstream>
#include <stdio.h>
#include <unistd.h>

namespace avro
{

namespace programs
{

using namespace graphics;

bool animating, a_pressed;
index_t nb_frames;
index_t nb_particles;
index_t nb_properties;
std::vector<gl_float> coordinates;
coord_t dim;
real_t dt;

gl_index vertex_array;

std::vector< std::string > macros = {"#version 410"};
std::shared_ptr<ShaderProgram> shader;

void
draw( Window& window ) {

  const mat4& view_matrix = window.camera().view_matrix();
  const mat4& projection_matrix = window.camera().projection_matrix();
  mat4 mvp = projection_matrix * view_matrix * window.plot(0).model_matrix();
  shader->setUniform("u_ModelViewProjectionMatrix" , mvp );

  if (animating) {

    GL_CALL( glBindVertexArray(vertex_array) );

    for (index_t i = 0; i < nb_frames; i++) {

      clock_t t0 = clock();

      // bind the data to the points buffer
      GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * nb_particles * 3 , coordinates.data() + 3*nb_particles*i , GL_STATIC_DRAW) );

      clock_t t1 = clock();
      real_t elapsed = real_t(t1 - t0)/real_t(CLOCKS_PER_SEC);

      // sleep for the desired time
      if (elapsed < dt) {
        usleep( 1e6 * (dt - elapsed) );
      }

      // clear the canvas
      glClearColor(1.0,1.0,1.0, 0.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // draw the fluid particles
      // if we don't use too much memory (maybe max 10k or 100k points with 1000 time steps ~ a few Gb), then maybe we can load everything into the same buffer
      GL_CALL( glDrawArrays( GL_POINTS , 0 , nb_particles ) );

      glfwSwapBuffers( window.window() );
    }
  }
  else {
    // clear the canvas
    glClearColor(1.0,1.0,1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // draw the fluid particles
    // if we don't use too much memory (maybe max 10k or 100k points with 1000 time steps ~ a few Gb), then maybe we can load everything into the same buffer
    GL_CALL( glDrawArrays( GL_POINTS , 0 , nb_particles ) );

    glfwSwapBuffers( window.window() );
  }

}

int
animate( int nb_input , const char** inputs )
{
#if AVRO_WITH_GL

  if (nb_input<1 || nb_input==-1)
  {
    printf("\t\tanimate [points file] [optional]\n");
    return 1;
  }

  const char **options = inputs +1;
  int nb_options = nb_input -1;

  ProcessCPU::initialize();
  ProcessGPU::initialize();

  // use float instead of double to represent coordinates
  using json_flt = nlohmann::basic_json<std::map, std::vector, std::string, bool,std::int64_t, std::uint64_t, float>;

  // read the json with the point data
  json_flt jin;
  std::ifstream file_in(inputs[0]);
  if (!file_in.good()) avro_assert_not_reached;
  file_in >> jin;

  // extract information about the points
  nb_frames = jin["nb_frames"];
  nb_particles = jin["nb_particles"];
  nb_properties = jin["nb_properties"];
  dim = jin["dimension"];

  printf("--> reading animation data: nb_frames = %lu, nb_particles = %lu, dim = %u\n",nb_frames,nb_particles,dim);

  std::vector<float> data = jin["coordinates"];
  std::vector<float> properties = jin["properties"];
  std::vector<gl_float> density = jin["density"];
  avro_assert( data.size() == nb_frames * nb_particles * dim );

  // process the coordinates
  index_t j = 0;
  coordinates.resize( nb_frames * nb_particles * 3 , 0.0 );
  for (index_t k = 0; k < nb_frames; k++) {
    for (index_t i = 0; i < nb_particles; i++)
    for (coord_t d = 0; d < dim; d++)
      coordinates[ k * (3*nb_particles) + i * 3 + d ] = gl_float(data[j++]);
  }

  // TODO process properties

  // determine the time step from the optional arguments
  bool found = false;
  dt = 0.1;
  found = parse<real_t>(lookfor(options,nb_options,"dt"),dt);
  printf("dt = %g\n",dt);

  // initialize the window
  int width = 600, height = width;
  Window window(width,height);
  window.init();

  // add a plot which won't be rendered, but we will use its model matrix
  CKF_Triangulation tri({2,2});
  Plot plot(tri);
  window.add_plot(&plot);

  window.compute_view();
  window.enable_controls(true);

  GL_CALL( glGenVertexArrays( 1, &vertex_array ) );
  GL_CALL( glBindVertexArray(vertex_array) );

  // create the points buffer
  gl_index coordinate_buffer;
  GL_CALL( glGenBuffers( 1 , &coordinate_buffer ) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , coordinate_buffer ) );

  // bind which attributes we want to draw
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
  GL_CALL( glEnableVertexAttribArray(0) );

  // create the density buffer
  gl_index density_buffer;
  GL_CALL( glGenBuffers( 1 , &density_buffer ) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , density_buffer ) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * nb_particles , density.data() , GL_STATIC_DRAW) );

  // bind which attributes we want to draw
  GL_CALL( glVertexAttribPointer( 1, 1, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 1 which holds density
  GL_CALL( glEnableVertexAttribArray(1) );

  shader = std::make_shared<ShaderProgram>("particles",false,macros);
  shader->use();

  glfwSetInputMode(window.window(), GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , coordinate_buffer ) );

  // bind the data for the first frame to the coordinates buffer
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, sizeof(gl_float) * nb_particles * 3 , coordinates.data() , GL_STATIC_DRAW) );

  vec3 center;
  float xmin = 1e20, xmax = -1e20;
  for (index_t k = 0; k < coordinates.size()/3; k++) {
    for (coord_t d = 0; d < 3; d++) {
      float xd = coordinates[3*k+d];
      center(d) += xd;
      if (xd < xmin) xmin = xd;
      if (xd > xmax) xmax = xd;
    }
  }

  for (coord_t d = 0; d < 3; d++)
    center(d) *= 3.0/coordinates.size();
  window.compute_view( center , 2.0*(xmax-xmin) );

  window.set_draw_callback( draw );

  // render loop
  while (true) {

    window.draw();
    glfwWaitEvents();

    // check if the user hit and released the A button, then we animate
    a_pressed = false;
    animating = false;
    if (glfwGetKey(window.window(), GLFW_KEY_A ) == GLFW_PRESS) a_pressed = true;
    if (a_pressed && glfwGetKey(window.window(), GLFW_KEY_A ) == GLFW_RELEASE) animating = true;

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window.window())) break;
    if (glfwGetKey(window.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }

#else
  printf("OpenGL is required to render particles\n");
#endif

  return 0;
}


} // programs

} // avro
