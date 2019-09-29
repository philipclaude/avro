#include "unit_tester.hpp"

#include "graphics/window.h"

using namespace ursa;
using namespace ursa::graphics;

UT_TEST_SUITE( WindowSuite )

UT_TEST_CASE( test1 )
{
  Window window( "ursa plotter" );

  window.initialize();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( WindowSuite )
