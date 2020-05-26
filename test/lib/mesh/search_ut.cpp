#include "unit_tester.hpp"

#include "library/ckf.h"

#include "mesh/search.h"

#include "numerics/geometry.h"

using namespace avro;

template<typename type>
static index_t
get_guess( const Topology<type>& topology , index_t skip )
{
  index_t guess = 0;
  while (topology.ghost(guess) || guess==skip)
    guess++;
  return guess;
}

UT_TEST_SUITE(mesh_search_suite)

UT_TEST_CASE(element_search_test)
{
  for (coord_t number=2;number<=4;number++)
  {
    coord_t dim = number;
    std::vector<index_t> dims(number,4);

    CKF_Triangulation topology(dims);
    topology.orient();
    topology.close();
    topology.neighbours().compute();
    topology.inverse().build();

    ElementSearch<Simplex> searcher(topology);

    // check the centroid of every element
    std::vector<real_t> ck(dim,0.);
    for (index_t k=0;k<topology.nb();k++)
    {
      if (topology.ghost(k)) continue;

      numerics::centroid( topology(k) , topology.nv(k) , topology.points() , ck );

      // find some other element to start from
      int ielem = searcher.find( ck.data() , get_guess(topology,k) );

      UT_ASSERT_EQUALS( ielem , int(k) );
    }

    // find a ghost and catch and exception
    index_t guess = 0;
    while (!topology.ghost(guess))
      guess++;
    UT_CATCH_EXCEPTION(searcher.find(ck.data(),guess));

    // loop through boundary facets
    std::vector<index_t> f(number,0);
    std::vector<real_t> cf(dim,0);
    std::vector<real_t> v(dim,0.);
    std::vector<real_t> p(dim,0.);
    std::vector<real_t> alpha(dim+1,0.);
    for (index_t k=0;k<topology.nb();k++)
    {
      if (!topology.ghost(k)) continue;

      avro_assert( topology(k,0)==0 );

      for (index_t j=0;j<number;j++)
        f[j] = topology(k,j+1);

      // get the centroid of this facet
      numerics::centroid( f.data() , f.size() , topology.points() , cf );

      // get the opposite neighbour
      int k1 = -1;
      for (index_t j=0;j<topology.neighbours().nfacets();j++)
      {
        if (topology.ghost(topology.neighbours()(k,j)))
          continue;
        k1 = topology.neighbours()(k,j);
        break;
      }
      UT_ASSERT( k1 >= 0 );

      // with exact predicates, this should be the element
      int ielem = searcher.find( cf.data() , get_guess(topology,k1) );
      UT_ASSERT_EQUALS( ielem , int(k1) );

      // get the index of k1
      int q = topology.neighbours().opposite(k1,k);
      avro_assert( q>=0 );

      // get a point outside the facet
      numerics::vector( topology.points()[q] , cf.data() , dim , v.data() );
      numerics::axpb( 1e-3 , v.data() , cf.data() , dim , p.data() );

      // make sure the search returns "outside"
      ielem = searcher.find( p.data() , get_guess(topology,k1) );
      UT_ASSERT_EQUALS( ielem , -1 );

      // even the brute-force search should return "outside"
      ielem = searcher.brute( p.data() );
      UT_ASSERT_EQUALS( ielem , -1 );

      // find the closest point and make sure it's the same simplex
      ielem = searcher.closest( p.data() , alpha );
      UT_ASSERT_EQUALS( ielem , k1 );
    }
  }

}
UT_TEST_CASE_END(element_search_test)

UT_TEST_SUITE_END(mesh_search_suite)
