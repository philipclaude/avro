#include "unit_tester.hpp"

#include "graphics/window.h"

using namespace avro;
using namespace avro::graphics;

UT_TEST_SUITE( WindowSuite )

UT_TEST_CASE( test1 )
{
  // must be run before the window starts
  glfwInit();

  Window window( "avro plot" , NULL );

  //window.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( WindowSuite )
