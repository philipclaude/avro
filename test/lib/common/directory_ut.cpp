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

#include "common/directory.h"
#include "avro_types.h"

using namespace avro;

UT_TEST_SUITE(directory_test_suite)

bool
contains( const std::vector<json>& listing , const std::string& d )
{
  for (index_t k=0;k<listing.size();k++)
    if (listing[k]["entry"]==d) return true;
  return false;
}

UT_TEST_CASE(test1)
{
  Directory dir(".");

  std::vector<json> listing;
  dir.ls( listing );

  //UT_ASSERT( contains(listing,"bin"));
  UT_ASSERT( contains(listing,"lib"));
  UT_ASSERT( contains(listing,"library"));
  UT_ASSERT( contains(listing,"regression"));
  //UT_ASSERT( contains(listing,"sandbox"));
  UT_ASSERT( contains(listing,"tmp"));

  UT_ASSERT( !contains(listing,"something-random") );

  dir.cd( "library" );
  listing.clear();
  dir.ls( listing );
  UT_ASSERT( contains(listing,"geometry"));
  UT_ASSERT( contains(listing,"meshes"));

  std::string ext = get_file_ext( dir.pwd() + "CMakeLists.txt" );
  UT_ASSERT_EQUALS( ext , "txt" );

}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(directory_test_suite)
