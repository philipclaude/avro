//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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

  std::shared_ptr<Object_t> prim;
  UT_CATCH_EXCEPTION( prim = std::make_shared<Object_t>(context,obj) );

  std::shared_ptr<EGADS::Body> body_ptr;
  body_ptr = std::make_shared<EGADS::Body>( context , obj );

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
