#include "unit_tester.hpp"

#include "geometrics/egads.h"
#include "geometrics/primitive.h"

#include "numerics/coordinate.h"

using namespace ursa;

UT_TEST_SUITE(PrimitiveSuite)

UT_TEST_CASE(test1)
{
  typedef geometrics::EGADS::Object GObject;

  GObject obj;
  geometrics::Primitive<GObject> prim(obj);

  numerics::Coordinate x(3);
  prim.project(x);

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(PrimitiveSuite)
