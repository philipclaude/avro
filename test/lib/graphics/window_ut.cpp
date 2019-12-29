#include "unit_tester.hpp"

#include "graphics/window.h"

using namespace luma;
using namespace luma::graphics;

UT_TEST_SUITE( WindowSuite )

UT_TEST_CASE( test1 )
{
  // must be run before the window starts
  glfwInit();

  Window window( "luma plot" , NULL );

  //window.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( WindowSuite )
