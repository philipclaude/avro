// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "common/mpi.hpp"
#include "common/parallel_for.h"
#include "common/process.h"
#include "common/tools.h"

using namespace avro;

#ifdef AVRO_MPI

void
main_hello( index_t k )
{
  printf("hello from main!\n");
}

void
worker_do( index_t k )
{
  printf("hello from worker %lu\n",k);
  avro_assert(k!=0);
}

void
main_end( index_t k )
{
  printf("goodbye!\n");
}

UT_TEST_SUITE(ParallelMPISuite)

UT_TEST_CASE(hello_test)
{
  ProcessMPI::main_begin( &main_hello );
  ProcessMPI::worker_do( &worker_do );
  //ProcessMPI::main_end( &main_end );
}
UT_TEST_CASE_END(hello_test)

class TestClass
{
public:
  typedef TestClass thisclass;

  TestClass() {}

  void main_hello( index_t k )
  {
    printf("hellooo from main\n");
    avro_assert( k==0 );
  }

  void worker_hello( index_t rank )
  {
    printf("hellooo from worker %lu\n",rank);

    if (ProcessMPI::nb_processes()>1) avro_assert( rank!=0 );
  }

  void main_goodbye( index_t rank )
  {
    printf("goodbye!\n");
  }

  void run()
  {
    ProcessMPI::main_begin( task_member( this , &thisclass::main_hello ) );
    ProcessMPI::worker_do( task_member( this , &thisclass::worker_hello ) );
    ProcessMPI::main_end( task_member( this , &thisclass::main_goodbye ) );
  }
};

UT_TEST_CASE(send_recv)
{
  TestClass T;
  T.run();

  real_t value = 1.;
  std::vector<int> values = {1,2,1,4,5};

  int tag1 = 0;
  int tag2 = 1;

  if (mpi::rank()==0)
  {

    for (index_t i=1;i<ProcessMPI::nb_processes();i++)
    {
      mpi::send( mpi::blocking{} , value , i  , tag1 );
      mpi::send( mpi::blocking{} , values , i , tag2 );
    }
  }
  else
  {
    real_t i = mpi::receive<real_t>( 0 , tag1 );
    UT_ASSERT_EQUALS( i , value );

    std::vector<int> x = mpi::receive< std::vector<int> >(0,tag2);
    UT_ASSERT_EQUALS( x.size() , values.size() );
    for (index_t j=0;j<values.size();j++)
    {
      printf("values[%lu] = %d\n",j,x[j]);
      UT_ASSERT_EQUALS( x[j] , values[j] );
    }
  }
}
UT_TEST_CASE_END(send_recv)

UT_TEST_CASE(bcast_test)
{
  UT_ASSERT_EQUALS( mpi::broadcast(10) , 10 );
  std::vector<real_t> x = {1.2,3.4};
  std::vector<real_t> y = mpi::broadcast(x);
  if (mpi::rank()!=0)
  {
    for (index_t i=0;i<x.size();i++)
      UT_ASSERT_NEAR( x[i] , y[i] , 1e-12 );
  }

  if (mpi::rank() == 0)
  {
    std::vector<int> data_to_send = {0,1,2,3,4,5,6,7,8,9};

    std::vector<mpi::request> requests;

    requests.emplace_back( mpi::send (mpi::async{}, data_to_send, 1) );

    // Computation happens here

    std::vector<mpi::status> const statuses = mpi::wait_all(requests);

    // Now std::vector<int> is deleted after the wait_all
  }
  else if (mpi::rank() == 1)
  {
    index_t j = 0;
    std::vector<int> x = mpi::receive<std::vector<int>>(0);
    for (index_t i=0;i<x.size();i++)
    {
      UT_ASSERT(i == j++);
    }
  }
}
UT_TEST_CASE_END(bcast_test)

UT_TEST_CASE(reduce_test)
{
  avro_assert( TEST_NUM_PROCS == 4 );

  UT_ASSERT(mpi::all_reduce(1, mpi::sum{}) == mpi::size());
  UT_ASSERT(mpi::all_reduce(1, mpi::prod{}) == 1);
  UT_ASSERT(mpi::all_reduce(mpi::rank(), mpi::max{}) == mpi::size() - 1);
  UT_ASSERT(mpi::all_reduce(mpi::rank(), mpi::min{}) == 0);

  std::vector<int> vec = {1,1,1};

  //int s = mpi::all_reduce(vec,mpi::sum{});
  for (auto const& sum : mpi::all_reduce(vec, mpi::sum{}))
  {
    UT_ASSERT(sum == mpi::size());
  }
  for (auto const prod : mpi::all_reduce(vec, mpi::prod{}))
  {
    UT_ASSERT(prod == 1);
  }

  UT_ASSERT(mpi::all_reduce(mpi::rank(), mpi::max{}) == mpi::size() - 1);
  UT_ASSERT(mpi::all_reduce(mpi::rank(), mpi::min{}) == 0);

  int root_sum = mpi::reduce(1, mpi::sum{}, 1);
  if (mpi::rank() == 1)
  {
    UT_ASSERT(root_sum == mpi::size());
  }

  auto const root_prod = mpi::reduce(1, mpi::prod{}, 1);
  if (mpi::rank() == 1)
  {
    UT_ASSERT(root_prod == 1);
  }

  auto const root_max = mpi::reduce(mpi::rank(), mpi::max{}, 1);
  if (mpi::rank() == 1)
  {
    UT_ASSERT(root_max == mpi::size() - 1);
  }

  auto const root_min = mpi::reduce(mpi::rank(), mpi::min{}, 1);
  if (mpi::rank() == 1)
  {
    UT_ASSERT(root_min == 0);
  }

  std::vector<int> array(2, mpi::rank());

  int host = 1;

  std::vector<int> root_sum2 = mpi::reduce(array, mpi::sum{}, host);
  print_inline( root_sum2 );
  if (mpi::rank() == host)
  {
    for (index_t i=0;i<root_sum2.size();i++)
    {
      UT_ASSERT_EQUALS(root_sum2[i],int(2*(ProcessMPI::nb_processes()-1)));
    }
  }

  std::vector<int> root_prod2 = mpi::reduce(array, mpi::prod{}, host);
  if (mpi::rank() == host )
  {
    for (index_t i=0;i<root_prod2.size();i++)
    {
      UT_ASSERT_EQUALS( root_prod2[i] , 0);
    }
  }

  std::vector<int> root_max2 = mpi::reduce(array, mpi::max{}, host);
  if (mpi::rank() == host)
  {
    for (index_t i=0;i<root_max2.size();i++)
    {
      UT_ASSERT_EQUALS(root_max2[i] , int(ProcessMPI::nb_processes()-1));
    }
  }

  std::vector<int> root_min2 = mpi::reduce(array, mpi::min{}, host);
  if (mpi::rank() == host)
  {
    for (index_t i=0;i<root_min2.size();i++)
    {
      UT_ASSERT_EQUALS( root_min2[i] , 0);
    }
  }
}
UT_TEST_CASE_END(reduce_test)

UT_TEST_SUITE_END(ParallelMPISuite)

#else // NO_MPIRUN

UT_TEST_SUITE(ParallelMPISuite)
UT_TEST_SUITE_END(ParallelMPISuite)

#endif // NO_MPIRUN
