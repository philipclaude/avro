#include "unit_tester.hpp"

#include "geometry/egads/body.h"
#include "geometry/egads/context.h"
#include "geometry/egads/object.h"
#include "geometry/entity.h"

#include "library/egads.h"

using namespace avro;

UT_TEST_SUITE(body_suite)

UT_TEST_CASE(test1)
{
  typedef EGADS::Object Object_t;
  EGADS::Context context;
  ego obj = nullptr;
  std::shared_ptr<Object_t> prim = std::make_shared<Object_t>(context,&obj);

  std::shared_ptr<EGADS::Body> body_ptr;
  UT_CATCH_EXCEPTION( body_ptr = std::make_shared<EGADS::Body>( context , prim->object() ) );

  //body.add( prim );
}
UT_TEST_CASE_END(test1)

UT_TEST_CASE(test2)
{
  EGADS::Context context;
  EGADS::Cube box(&context,{1,1,1});

  //box.add( prim );
}
UT_TEST_CASE_END(test2)

UT_TEST_SUITE_END(body_suite)
