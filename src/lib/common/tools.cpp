//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/mpi.hpp"
#include "common/process.h"
#include "common/tools.h"

#include "graphics/gl.h"

#include "numerics/predicates.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <math.h>
#include <cstdlib>

namespace avro
{

real_t
random_within( const real_t lo , const real_t hi )
{
	real_t s = (real_t) rand() / (real_t) RAND_MAX;
  return lo +s*(hi-lo);
}

int
random_within( const int lo , const int hi )
{
	real_t s = (real_t)rand()/(real_t) RAND_MAX;
	return lo +s*(hi -lo);
}


void
initialize_avro()
{
  printf("\n================================================\n");
  printf(  "| avro -- (c) Philip Claude Caplan (2017-2021) |\n");
  printf("================================================\n");

  // initialize the process managers
  ProcessMPI::initialize();
  ProcessCPU::initialize();
  ProcessGPU::initialize();
  #if AVRO_MPI
  printf("--> initialized MPI manager: %s (mpi::size = %lu).\n",ProcessMPI::manager_name().c_str(),mpi::size());
  #endif
  printf("--> initialized CPU manager: %s (nb_thread = %lu).\n",ProcessCPU::manager_name().c_str(),ProcessCPU::maximum_concurrent_threads());
  printf("--> initialized GPU manager: %s (nb_cores  = %lu).\n",ProcessGPU::manager_name().c_str(),ProcessGPU::maximum_concurrent_threads());


  // initialize the predicates
  exactinit(0,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();
  printf("--> initialized robust predicates.\n");

	#if 0
  if (!glfwInit()) printf("error initializing glfw!\n");

  // set the version
  #if AVRO_HEADLESS_GRAPHICS // core 3.3 supported by wazowski's drivers
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
  #else
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  #endif
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow* window = glfwCreateWindow( 1,1,"avro" , NULL, NULL);
  avro_assert( window!=NULL );
  glfwMakeContextCurrent(window);
  gladLoadGL();

  const GLubyte *renderer = glGetString( GL_RENDERER );
  const GLubyte *glslVersion = glGetString( GL_SHADING_LANGUAGE_VERSION );
  GLint major,minor;
  glGetIntegerv(GL_MAJOR_VERSION, &major);
  glGetIntegerv(GL_MINOR_VERSION, &minor);
  printf("--> initialized OpenGL: %s (version %d.%d with GLSL %s)\n",renderer,major,minor,glslVersion);
	#endif

}

}
