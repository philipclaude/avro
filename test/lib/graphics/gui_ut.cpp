//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "graphics/gl.h"
#include "graphics/gui.h"
#include "graphics/window.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( gui_test_suite )


UT_TEST_CASE( imgui_test )
{

  if (AVRO_FULL_UNIT_TEST) return;

  Window window(600,600);
  window.init();
  //window.set_enable_controls(false);

  GUI gui(window);

  window.draw();

  while (true) {

    // wait for user input
    glfwWaitEvents();

    gui.draw();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window.window())) break;
    if (glfwGetKey(window.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}
UT_TEST_CASE_END( imgui_test )

UT_TEST_SUITE_END( gui_test_suite )
