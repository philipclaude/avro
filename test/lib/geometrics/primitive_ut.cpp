#include "unit_tester.hpp"

#include "geometry/egads/context.h"
#include "geometry/egads/object.h"
#include "geometry/entity.h"

#include "numerics/coordinate.h"

using namespace luna;

UT_TEST_SUITE(PrimitiveSuite)

UT_TEST_CASE(test1)
{
  EGADS::Context context;
  ego obj;
  EGADS::Object prim(context,&obj);

  numerics::Coordinate x(3),u(2);
  prim.inverse(x,u);

  prim.evaluate(u,x);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(PrimitiveSuite)
