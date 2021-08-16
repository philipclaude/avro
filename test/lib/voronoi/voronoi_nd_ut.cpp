#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/optimal_transport.h"

#include <json/json.hpp>

#include <fstream>

using namespace avro;

UT_TEST_SUITE( voronoi_cell_test_suite )

UT_TEST_CASE( test_nd )
{
  coord_t numberL = 2;
  coord_t numberH = 4;
  coord_t powerL = 2;
  coord_t powerH = 2;
  bool generate_mode = false;

  for (coord_t number = numberL; number <= numberH; number++)
  {
    coord_t dim = number;

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

    //ckf.points().incidence().print();

    // create the polytopal topology
    Topology<Polytope> domain(ckf.points(),number);

    std::vector<index_t> cube = linspace( std::pow(2,dim) );
    print_inline(cube);
    domain.add( cube.data() , cube.size() );
    for (coord_t power = powerL; power <= powerH; power++)
    {
      index_t nb_points = std::pow( 2 , power );
      printf("dim = %u, nb_points = %lu\n",dim,nb_points);

      // create random delaunay vertices
      Points delaunay( dim );
      std::vector<real_t> p(dim,0.);
      for (index_t k=0;k<nb_points;k++)
      {
        for (index_t d=0;d<dim;d++)
          p[d] = random_within(0.,1.);
        delaunay.create(p.data());
      }

      index_t nb_iter = 20;
      for (index_t iter = 0; iter <= nb_iter; iter++)
      {
        voronoi::LaguerreDiagram<Polytope> diagram( delaunay , domain );
        diagram.compute(false);

        for (index_t k=0;k<diagram.nb();k++)
        {
          real_t dmin =  1e20;
          real_t dmax = -1e20;

          // approximate the centroid
          std::vector<real_t> c( dim , 0.0 );
          for (index_t j=0;j<diagram.nv(k);j++)
          {
            real_t r = numerics::distance( diagram.points()[ diagram(k,j) ] , delaunay[k] , dim );
            if (r < dmin) r = dmin;
            if (r > dmax) r = dmax;
            for (coord_t d=0;d<dim;d++)
              c[d] += r*diagram.points()[ diagram(k,j) ][d];
          }

          for (coord_t d=0;d<dim;d++)
            delaunay[k][d] = (c[d] / diagram.nv(k)) / dmax;
        }

        if (iter == nb_iter)
        {

          // save to a file
          #if 0
          json J;
          J["dim"] = dim;
          J["nb"]  = delaunay.nb();

          std::vector<real_t> x( delaunay.nb() * delaunay.dim() );
          index_t i = 0;
          for (index_t k=0;k<delaunay.nb();k++)
          for (index_t d=0;d<delaunay.dim();d++)
            x[i++] = delaunay[k][d];
          J["x"] = x;

          std::ofstream output("samples-dim"+std::to_string(dim)+"-n"+std::to_string(delaunay.nb())+".json");
          output << std::setw(4) << J << std::endl;
          #endif

          /*
          if (number > 3 || (nb_points >= 1e4) || generate_mode) continue;
          graphics::Viewer vis;
          vis.add(diagram);
          vis.run(AVRO_FULL_UNIT_TEST);
          */
        }
      } // loop over power
    }
  } // loop over number

}
UT_TEST_CASE_END( test_nd )

UT_TEST_SUITE_END( voronoi_cell_test_suite )
