#include "unit_tester.hpp"

#include "library/obj.h"

using namespace avro;
using namespace avro::library;

UT_TEST_SUITE( obj_file_suite )

UT_TEST_CASE( test1 )
{

  objFile topology( BASE_TEST_DIR+"/geometry/obj/suzanne.obj" );

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( obj_file_suite )
