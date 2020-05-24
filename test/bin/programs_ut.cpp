#include "unit_tester.hpp"

#include "../bin/programs.h"

#include <string>
#include <vector>

using namespace avro;

UT_TEST_SUITE(avro_programs_suite)


UT_TEST_CASE(test1)
{
  int result;

  const char* command1[] = {"CKF-3-3-3","box","Linear-3d","tmp/cl","nb_iter=1"};
  result = programs::adapt(5,command1);
  UT_ASSERT_EQUALS( result , 0 );

  const char* command2[] = {"CKF-3-3-3","random","none","nb_iter=5"};
  result = programs::voronoi(4,command2);
  UT_ASSERT_EQUALS( result , 0 );

  const char* command3[] = {"CKF-3-3-3","sample","none","nb_iter=5"};
  result = programs::voronoi(4,command3);
  UT_ASSERT_EQUALS( result , 0 );

  const char* command4[] = {"CKF-3-3-3","Linear-3d","tmp/output.json","nb_expected=10"};
  result = programs::conformity(4,command4);
  UT_ASSERT_EQUALS( result , 1 );
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(avro_programs_suite)
