#include "unit_tester.hpp"

#include "types.h"

#include "blossom5/PerfectMatching.h"

using namespace avro;

UT_TEST_SUITE( blossom_test_suite )

UT_TEST_CASE( test1 )
{

  PerfectMatching graph( 16 , 53 );

  int infty = PM_INFTY/2;
  printf("PM infty = %d\n",PM_INFTY);

  graph.AddEdge( 0,3,321 );
  graph.AddEdge( 0,4,341 );
  graph.AddEdge( 0,9,338 );
  graph.AddEdge( 0,11,326 );
  graph.AddEdge( 0,12,infty );
  graph.AddEdge( 0,13,337 );
  graph.AddEdge( 0,15,347 );
  graph.AddEdge( 1,2,319 );
  graph.AddEdge( 1,5,318 );
  graph.AddEdge( 1,8,infty );
  graph.AddEdge( 1,10,318 );
  graph.AddEdge( 1,13,330 );
  graph.AddEdge( 1,14,324 );
  graph.AddEdge( 1,15,340 );
  graph.AddEdge( 2,4,333 );
  graph.AddEdge( 2,5,317 );
  graph.AddEdge( 2,6,324 );
  graph.AddEdge( 2,7,322 );
  graph.AddEdge( 2,8,321 );
  graph.AddEdge( 2,10,317 );
  graph.AddEdge( 2,11,infty );
  graph.AddEdge( 2,12,326 );
  graph.AddEdge( 2,13,329 );
  graph.AddEdge( 3,4,328 );
  graph.AddEdge( 3,5,infty );
  graph.AddEdge( 3,8,316 );
  graph.AddEdge( 3,10,312 );
  graph.AddEdge( 3,12,321 );
  graph.AddEdge( 4,7,337 );
  graph.AddEdge( 5,8,infty );
  graph.AddEdge( 5,13,328 );
  graph.AddEdge( 5,15,338 );
  graph.AddEdge( 6,7,infty );
  graph.AddEdge( 7,10,infty );
  graph.AddEdge( 7,12,330 );
  graph.AddEdge( 8,10,320 );
  graph.AddEdge( 8,11,321 );
  graph.AddEdge( 8,12,329 );
  graph.AddEdge( 8,13,332 );
  graph.AddEdge( 8,14,326 );
  graph.AddEdge( 9,11,330 );
  graph.AddEdge( 9,12,infty );
  graph.AddEdge( 9,14,infty );
  graph.AddEdge( 9,15,infty );
  graph.AddEdge( 10,12,325 );
  graph.AddEdge( 11,12,326 );
  graph.AddEdge( 11,13,329 );
  graph.AddEdge( 11,15,339 );
  graph.AddEdge( 12,13,337 );
  graph.AddEdge( 12,14,infty );
  graph.AddEdge( 12,15,347 );
  graph.AddEdge( 13,15,infty );
  graph.AddEdge( 14,15,344 );

  graph.Solve();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( blossom_test_suite )
