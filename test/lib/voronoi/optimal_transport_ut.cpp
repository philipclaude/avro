#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "voronoi/optimal_transport.h"

#include <fstream>

using namespace avro;

UT_TEST_SUITE( optimal_transport_test_suite )

UT_TEST_CASE( test_nd_polytope )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t numberL = 3;
  coord_t numberH = 3;

  index_t nb_points = 1e2;

  for (coord_t number = numberL; number <= numberH; number++) {

    coord_t dim = number;

    CubeDomain<type> domain(dim,dim,10);
    printf("dim = %u, nb_points = %lu\n",dim,nb_points);

    // create random delaunay vertices
    Points delaunay( dim );
    std::vector<index_t> elems;
    for (index_t k = 0; k < nb_points; k++) {

      index_t elem;
      std::vector<real_t> p(dim,0.);
      if (typeid(type) == typeid(Simplex)) {

        // pick a random element
        elem = random_within(int(0),int(domain.nb()));
        index_t N = domain.nv(elem);
        std::vector<real_t> alpha(N,0.0);

        // randomize the barycentric coordinates
        alpha[N-1] = 1.0;
        for (index_t j = 0; j < N-1; j++) {
          alpha[j]   = random_within(0.,1.0);
          alpha[N-1] -= alpha[j];
        }
        for (index_t j = 0; j < N; j++) {
          if (alpha[j] < 0.0) alpha[j] = 0.0;
        }

        // compute the point coordinates
        for (coord_t d = 0; d < dim; d++)
        for (index_t j = 0; j < N; j++)
          p[d] += alpha[j]*domain.points()[ domain(elem,j) ][d];
      }
      else {
        elem = 0;
        for (coord_t d = 0; d < dim; d++)
          p[d] = random_within(0.0,1.0);
      }

      // create the delaunay point and retain which element in the domain it is in
      delaunay.create(p.data());
      elems.push_back( elem );
    }

    // initialize and compute the laguerre diagram
    voronoi::LaguerreDiagram<type> diagram( delaunay , domain );
    diagram.set_elements( elems );

    voronoi::IntegrationSimplices triangulation(number,number);

    diagram.compute(false,&triangulation);

    std::shared_ptr<voronoi::TriangulationCells> tc = std::make_shared<voronoi::TriangulationCells>(triangulation);
    triangulation.fields().make("c",tc);

    std::shared_ptr<voronoi::TriangulationElements> te = std::make_shared<voronoi::TriangulationElements>(triangulation);
    triangulation.fields().make("e",te);

    #if 0
    if (number > 3 || (nb_points >= 1e6)) continue;
    graphics::Viewer vis;
    vis.add(triangulation);
    //vis.add(diagram);
    vis.run(AVRO_FULL_UNIT_TEST);
    #elif 1
    graphics::OpenGL_Application app;
    app.add( triangulation );
    app.run(AVRO_FULL_UNIT_TEST);
    #endif

  } // loop over number

}
UT_TEST_CASE_END( test_nd_polytope )

UT_TEST_SUITE_END( optimal_transport_test_suite )
