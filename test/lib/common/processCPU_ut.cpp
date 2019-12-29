#include "unit_tester.hpp"

#include "common/parallel_for.h"

using namespace luma;

class Threadable
{
public:
  typedef Threadable thisclass;

  Threadable() {}

  void hello_from_thread( index_t i )
  {
    printf("hello from thread %lu\n",i);
  }

  void run_in_parallel( const index_t nthread )
  {
    ProcessCPU::parallel_for (
      parallel_for_member_callback( this , &thisclass::hello_from_thread ),
      0,nthread
    );
  }
};

UT_TEST_SUITE(ParallelCPUSuite)

UT_TEST_CASE(threadable_test)
{
  // openmp thread manager
  Threadable T;

  T.run_in_parallel(10);
}
UT_TEST_CASE_END(threadable_test)

UT_TEST_SUITE_END(ParallelCPUSuite)
