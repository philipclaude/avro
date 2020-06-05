//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "graphics/application.h"
#include "graphics/shader.h"

#include "mesh/points.h"
#include "mesh/topology.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( graphics_shader_test_suite )

UT_TEST_CASE( uniforms )
{
  Visualizer vis;

  ShaderProgram& shader = vis.manager().shaders()["wv"];

  int loc = shader.getUniformLocation("test");
  UT_ASSERT_EQUALS( loc , -1 );

  loc = shader.getUniformLocation("lightDir");
  UT_ASSERT( loc >= 0 );
  shader.setUniform( "lightDir" , 1. , 1. , 1. );

  float v[3] = {1.,1.,1.};
  shader.setUniform( "lightDir" , 1 , v );

  UT_CATCH_EXCEPTION( shader.setUniform( "lightDir" , 1. , 1. , 1. , 1. ) );

  vec2 v2 = {1,1};
  UT_CATCH_EXCEPTION( shader.setUniform( "lightDir" , v2  ) );

  vec3 v3 = {1,1,1};
  shader.setUniform( "lightDir" , v3  );

  vec4 v4 = {1,1,1,1};
  UT_CATCH_EXCEPTION( shader.setUniform( "lightDir" , v4  ) );

  mat3 m3;
  UT_CATCH_EXCEPTION( shader.setUniform( "MVP" , m3  ) );

  mat4 mvp;
  shader.setUniform( "MVP" , mvp );

  shader.setUniform( "test" , true );


  shader.printActiveUniforms();
  shader.printActiveAttribs();
  shader.validate();

}
UT_TEST_CASE_END( uniforms )

UT_TEST_SUITE_END( graphics_shader_test_suite )
