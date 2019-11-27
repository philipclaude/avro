#include "unit_tester.hpp"

#include "library/obj.h"

using namespace luna;
using namespace luna::library;

UT_TEST_SUITE( objFile_suite )

UT_TEST_CASE( test1 )
{

  objFile topology( "/Users/pcaplan/Desktop/suzanne.obj" );

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( objFile_suite )
