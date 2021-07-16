#include "unit_tester.hpp"

#include "graphics/gl.h"
#include "graphics/new/gui.h"
#include "graphics/new/window.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( gui_test_suite )


UT_TEST_CASE( imgui_test )
{

  Window window(600,600);
  window.init();
  //window.set_enable_controls(false);

  GUI gui(window);
  window.set_gui(&gui);

  window.draw();

  while (true) {

    // wait for user input
    glfwWaitEvents();

    gui.begin_draw();
    gui.draw();

    // determine if we should exit the render loop
    if (glfwWindowShouldClose(window.window())) break;
    if (glfwGetKey(window.window(), GLFW_KEY_ESCAPE ) == GLFW_PRESS) break;
  }
}
UT_TEST_CASE_END( imgui_test )

UT_TEST_SUITE_END( gui_test_suite )
