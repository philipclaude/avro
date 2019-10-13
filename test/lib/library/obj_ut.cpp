#include "unit_tester.hpp"

#include "library/obj.h"

using namespace ursa;
using namespace ursa::library;

UT_TEST_SUITE( objFile_suite )

UT_TEST_CASE( test1 )
{

  objFile topology( "/Users/pcaplan/Desktop/suzanne.obj" );

  topology.printData();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( objFile_suite )