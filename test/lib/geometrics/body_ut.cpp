#include "unit_tester.hpp"

#include "geometrics/body.h"
#include "geometrics/egads.h"
#include "geometrics/primitive.h"

using namespace ursa;

UT_TEST_SUITE(BodySuite)

UT_TEST_CASE(test1)
{
  typedef geometrics::EGADS::Object Object_t;

  geometrics::EGADS::Context context;
  ego obj;
  std::shared_ptr<Object_t> prim = std::make_shared<Object_t>(context,&obj);

  geometrics::EGADS::Body body( context , prim->object() );

  body.add( prim );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(BodySuite)
