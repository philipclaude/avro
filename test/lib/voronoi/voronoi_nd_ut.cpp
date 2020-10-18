#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi_cell.h"

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_nd )
{
  coord_t number = 3;
  coord_t dim = number;
  index_t nb_points = std::pow( 2 , 10 );
  printf("nb_points = %lu\n",nb_points);

  // create a CKF triangulation with only 2 points in each direction
  std::vector<index_t> dims(number,2);
  CKF_Triangulation ckf(dims);

  // compute the facets of the mesh
  Facets facets(ckf);
  facets.compute();

  // determine which facets are on a common hyperplane
  std::vector<std::vector<int>> v2b( ckf.points().nb() );
  std::vector<index_t> f(number);
  std::vector<real_t*> x(number);
  std::vector<real_t> n(dim);
  std::vector<real_t> c(dim);
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,f);

    // retrieve the facet points
    for (index_t j=0;j<number;j++)
      x[j] = ckf.points()[ f[j] ];

    // compute the normal to the facet
    numerics::normal( x , n.data() , dim );

    // which direction is this in?
    coord_t dir = 0;
    for (coord_t d=0;d<dim;d++)
    {
      if (fabs( fabs(n[d])-1.0) < 1e-8)
      {
        dir = d;
        break;
      }
    }

    // are we at 0 or 1?
    numerics::centroid( f.data() , f.size() , ckf.points() , c );
    bool zero = true;
    if (fabs(c[dir]) > 1e-8)
    {
      avro_assert( fabs(c[dir] -1.0) < 1e-8 );
      zero = false;
    }

    // determine which plane this corresponds to
    int b;
    if (zero) b = - dir - 1;
    else b = - dim - dir - 1;

    // let each vertex know it is on this bisector
    for (coord_t d=0;d<f.size();d++)
    {
      v2b[f[d]].push_back( b );
    }
  }

  // set the incidence relations
  for (index_t k=0;k<ckf.points().nb();k++)
  {
    uniquify( v2b[k] );
    ckf.points().incidence().add( v2b[k].data() , v2b[k].size() );
  }

  ckf.points().incidence().print();

  // create the polytopal topology
  Topology<Polytope> domain(ckf.points(),number);

  std::vector<index_t> cube = linspace( std::pow(2,dim) );
  print_inline(cube);
  domain.add( cube.data() , cube.size() );

  // create random delaunay vertices
  Delaunay delaunay( dim );
  std::vector<real_t> p(dim,0.);
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t d=0;d<dim;d++)
      p[d] = random_within(0.,1.);
    delaunay.create(p.data());
  }

  delaunay::VoronoiDiagram diagram( delaunay , domain );
  diagram.compute(false);

  // option to visualize the voronoi diagram
  if (number > 3 || nb_points >= 1e4) return;
  graphics::Visualizer vis;
  vis.add_topology(diagram);
  vis.run();

}
UT_TEST_CASE_END( test_nd )

UT_TEST_SUITE_END( voronoi_cell_test_suite )
