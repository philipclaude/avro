#include "unit_tester.hpp"

#include "geometrics/egads.h"
#include "geometrics/primitive.h"

#include "numerics/coordinate.h"

using namespace ursa;

UT_TEST_SUITE(PrimitiveSuite)

UT_TEST_CASE(test1)
{
  geometrics::EGADS::Context context;
  ego obj;
  geometrics::EGADS::Object prim(context,&obj);

  numerics::Coordinate x(3),u(2);
  prim.inverse(x,u);

  prim.evaluate(u,x);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(PrimitiveSuite)
