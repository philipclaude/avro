#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include "measures.h"
#include "visualize.h"

#include <fstream>
#include <iomanip>

#include <sys/time.h>
#include <sys/resource.h>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

UT_TEST_CASE( test1 )
{
  typedef Polytope type;

    for (coord_t number = 2; number <= 4; number++)
    {
      std::vector<real_t> time_voronoi,time_neighbours,time_triangulation;
      std::vector<index_t> number_of_points;
      std::vector<real_t> time_integration[3];

      std::vector<index_t> nb_points(3);
      if (number == 2) {
        nb_points[0] = 1e4;
        nb_points[1] = 1e5;
        nb_points[2] = 1e6;
      }
      if (number == 3) {
        nb_points[0] = 1e4;
        nb_points[1] = 1e5;
        nb_points[2] = 2.5e5;
      }
      if (number == 4) {
        nb_points[0] = 1e3;
        nb_points[1] = 1e4;
        nb_points[2] = 1.5e4;
      }

      for (real_t i  = 0; i < nb_points.size(); i++)
      {
        coord_t dim = number+1;
        CubeDomain<type> domain(number,dim,2);

        index_t nb_point = nb_points[i];
        printf("--> dim = %u, %lu points\n",number,nb_point);

        delaunay::DensityMeasure_Uniform density(1.0);
        delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density);
        transport.sample( nb_point );
        transport.set_mode(0); // points mode

        transport.compute_laguerre();

        for (coord_t quad_order = 2; quad_order <= 4; quad_order++)
        {
          std::vector<real_t> dc_dx( nb_point*number );

          transport.quad_order() = quad_order;

          transport.evaluate( dc_dx.data() , nullptr );

          time_integration[quad_order-2].push_back( transport.time_integrate() );
        }

        number_of_points.push_back(nb_point);
        time_voronoi.push_back( transport.diagram().time_voronoi() );
        time_neighbours.push_back( transport.diagram().time_neighbours() );
        time_triangulation.push_back( transport.diagram().time_decompose() );
      }

      std::string filename_output = "/Users/pcaplan/Dropbox/research/publications/imr-2021-xxxx/performance/whitenoise-timing-breakdown-dim-" + std::to_string(number) + ".json";

      json J;
      J["nb_points"] = number_of_points;
      J["time_voronoi"] = time_voronoi;
      J["time_neighbours"] = time_neighbours;
      J["time_triangulation"] = time_triangulation;
      J["time_integration-2"] = time_integration[0];
      J["time_integration-3"] = time_integration[1];
      J["time_integration-4"] = time_integration[2];

      std::ofstream output(filename_output);
      output << std::setw(4) << J << std::endl;
      output.close();
  }

  // if (number > 3) return;
  //
  // graphics::Visualizer vis;
  // vis.add_topology( transport.diagram() );
  // vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
