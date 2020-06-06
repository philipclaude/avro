#include "unit_tester.hpp"

#include <stdio.h>

UT_TEST_SUITE( sandbox_philip_test1 )

UT_TEST_CASE( test1 )
{
  printf("hello!\n");
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_philip_test1 )
