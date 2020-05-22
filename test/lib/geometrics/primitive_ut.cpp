#include "unit_tester.hpp"

#include "geometry/egads/context.h"
#include "geometry/egads/object.h"
#include "geometry/entity.h"

#include "numerics/coordinate.h"

using namespace avro;

UT_TEST_SUITE( geometry_primitive_suite)

UT_TEST_CASE(test1)
{
  EGADS::Context context;
  ego obj = nullptr;
  std::shared_ptr<EGADS::Object> prim;
  UT_CATCH_EXCEPTION( prim = std::make_shared<EGADS::Object>(context,&obj) );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(geometry_primitive_suite)
