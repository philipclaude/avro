#include "unit_tester.hpp"

#include "geometry/egads/body.h"
#include "geometry/egads/context.h"
#include "geometry/egads/object.h"
#include "geometry/entity.h"

using namespace luna;

UT_TEST_SUITE(BodySuite)

UT_TEST_CASE(test1)
{
  typedef EGADS::Object Object_t;

  EGADS::Context context;
  ego obj;
  std::shared_ptr<Object_t> prim = std::make_shared<Object_t>(context,&obj);

  EGADS::Body body( context , prim->object() );

  body.add( prim );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(BodySuite)
