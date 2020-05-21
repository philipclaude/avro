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
  EGADS::Object prim(context,&obj);

  numerics::Coordinate x(3),u(2);
  prim.inverse(x,u);

  prim.evaluate(u,x);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(geometry_primitive_suite)
