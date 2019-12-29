#include "unit_tester.hpp"

#include "library/ckf.h"
#include "library/metric.h"

#include "adaptation/filter.h"
#include "adaptation/metric.h"
#include "adaptation/primitive.h"

using namespace luma;

UT_TEST_SUITE( insertion_filter_suite )

UT_TEST_CASE(filter3d)
{
  const coord_t number = 3;
  std::vector<index_t> dims( number , 10 );
  CKF_Triangulation cube( dims );

  const coord_t dim = cube.points().dim();

  // add the cube points as fixed points
  Filter filter( cube.points().dim() );
  for (index_t k=0;k<cube.points().nb();k++)
    filter.createPermanent( cube.points()[k] );

  // add some points along the edges
  std::vector<index_t> edges;
  cube.get_edges(edges);
  for (index_t k=0;k<edges.size()/2;k++)
  {
    real_t s = random_within( 0. , 1. );
    std::vector<real_t> x(dim);
    index_t p0 = edges[2*k];
    index_t p1 = edges[2*k+1];
    real_t* x0 = cube.points()[p0];
    real_t* x1 = cube.points()[p1];
    for (coord_t d=0;d<dim;d++)
      x[d] = x0[d] +s*(x1[d]-x0[d]);

    std::vector<real_t> u(number-1,0);
    filter.createCandidate( p0 , p1 , s , x.data() , NULL , u.data() );
  }
  UT_ASSERT_EQUALS( filter.nb() , cube.points().nb()+edges.size()/2 );
  UT_ASSERT_EQUALS( filter.nb_candidates() , edges.size()/2 );

  // part 2
  Insert<Simplex> inserter(cube);
  library::MetricField_Uniform analytic(dim,2);
  MetricAttachment attachment(analytic,cube.points());
  MetricField<Simplex> metric(cube,attachment);
  filter.generateCandidates( cube , metric , 2. , inserter );
  filter.buildTree();

  printf("nb insertion candidates = %lu\n",filter.nb_candidates());

  // add the candidates
  index_t idx = cube.points().nb();
  std::vector<index_t> insertions;
  for (index_t k=0;k<filter.nb_candidates();k++)
  {

    // remember to insert this vertex for later
    insertions.push_back( filter.candidate(k) );

    // accept the candidate
    filter.accept( filter.candidate(k) ,idx++ );
  }
}
UT_TEST_CASE_END(filter3d)

UT_TEST_SUITE_END(insertion_filter_suite)
