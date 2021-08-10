#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/new/cell.h"

#include <nnsearch/nn_search.h>

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_2d )
{
  GEO::PCK::initialize();

  static coord_t number = 3;
  static coord_t dim = number;
  index_t nb_points = 1e1;

  std::vector<index_t> dims(number,2);
  CKF_Triangulation ckf(dims);

  Points points(dim);
  real_t x0[] = {0,0,0}; points.create(x0);
  real_t x1[] = {1,0,0}; points.create(x1);
  real_t x2[] = {0,1,0}; points.create(x2);
  real_t x3[] = {0,0,1}; points.create(x3);

  Topology<Simplex> topology(points,number);
  index_t t[] = {0,1,2,3};
  topology.add(t,number+1);
  topology.build_structures();

  Points sites(dim);
  ckf.points().copy( sites );
  nb_points = 1;//sites.nb();

  // all points are clipped against element 0
  index_t elem = 0;

  GEO::NearestNeighborSearch* searcher = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> X(sites.nb()*dim);
  for (index_t k = 0; k < sites.nb(); k++) {
    for (index_t d = 0; d < dim; d++)
      X[k*dim+d] = sites[k][d];
  }
  searcher->set_points( sites.nb() , X.data() );

  Points vertices(dim);
  Topology<Polytope> diagram(vertices,number);

  printf("computing cells\n");
  for (index_t k = 0; k < nb_points; k++) {

    printf("computing voronoi cell for site %lu\n",k);
    sites.print(k);

    voronoi::Cell cell( k , sites , topology , *searcher );
    cell.compute(elem);
    cell.Table<index_t>::print();
    cell.points().print();

    index_t nb_points = vertices.nb();
    for (index_t j = 0; j < cell.points().nb(); j++) {
      vertices.create( cell.points()[j] );
      std::vector<int> b = cell.points().incidence().get(j);
      vertices.incidence().add( b.data() , b.size() );
    }

    for (index_t j = 0; j < cell.nb(); j++) {
      std::vector<index_t> polytope = cell.get(j);
      for (index_t i = 0; i < polytope.size(); i++)
        polytope[i] += nb_points;
      diagram.add( polytope.data() , polytope.size() );
    }

//    vis.add(cell);
  }
  printf("done computing cells\n");

  graphics::Viewer vis;
  vis.add(diagram);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test_2d )


UT_TEST_SUITE_END( voronoi_cell_test_suite )
