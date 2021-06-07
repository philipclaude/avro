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

  for (index_t inoise = 0; inoise < 2; inoise++)
  {
    bool white = false;
    if (inoise == 0) white = true;
    else white = false;

    for (coord_t number = 2; number <= 6; number++)
    {
      std::vector<real_t> time_voronoi,time_neighbours,time_triangulation;
      std::vector<index_t> unshared_data,unshared_stack;
      std::vector<index_t> number_of_points;
      std::vector<real_t> number_of_vertices,number_of_facets;

      coord_t pmax = 6;
      if (number > 4) pmax = 5;
      if (number > 5) pmax = 4;

      for (coord_t power = 3; power <= pmax; power++)
      {

        std::string filename = "/home/pcaplan/Dropbox/research/publications/imr-2021-xxxx/performance/points-dim-" + std::to_string(number) + "-p-" + std::to_string(power) + ".json";
        printf("filename = %s\n",filename.c_str());
        std::fstream file;
        file.open(filename);
        avro_assert( file.is_open() );
        std::stringstream s;
        s << file.rdbuf();

        json J = json::parse( s.str().c_str() );

        std::vector<real_t> coordinates = J.at("points");

        coord_t dim = number+1;
        CubeDomain<type> domain(number,dim,2);

        index_t nb_points = coordinates.size() / number;
        printf("--> dim = %u, imported %lu points\n",number,nb_points);

        delaunay::SemiDiscreteOptimalTransport<type> transport(domain,nullptr);
        transport.sample( nb_points );

        // if blue noise, use the actual samples
        if (!white)
          transport.set_delaunay( coordinates.data() , number );

        transport.compute_laguerre();

        const LaguerreDiagram<type>& laguerre = transport.diagram();

        index_t nbv = 0, nbf = 0;
        std::vector<int> hrep;
        for (index_t k = 0; k < laguerre.nb(); k++)
        {
          // assuming 1 cell = 1 power cell
          nbv += laguerre.nv(k);

          laguerre.element().hrep( laguerre(k) , laguerre.nv(k) , hrep );
          nbf += hrep.size();
        }


        // this doesn't seem to work :/
        rusage usage;
        getrusage( RUSAGE_SELF , &usage );
        //printf("usage = %lu %lu\n",usage.ru_idrss,usage.ru_isrss);
        unshared_data.push_back( usage.ru_idrss );
        unshared_stack.push_back( usage.ru_isrss );


        number_of_points.push_back(nb_points);
        time_voronoi.push_back( transport.diagram().time_voronoi() );
        time_neighbours.push_back( transport.diagram().time_neighbours() );
        time_triangulation.push_back( transport.diagram().time_decompose() );
        number_of_vertices.push_back( real_t(nbv) / laguerre.nb() );
        number_of_facets.push_back( real_t(nbf) / laguerre.nb() );
      }

      std::string noise_type = ((white) ? "white" : "blue");

      std::string filename_output = "/home/pcaplan/Dropbox/research/publications/imr-2021-xxxx/performance/timing-" + noise_type + "-dim-" + std::to_string(number) + ".json";

      json J;
      J["nb_points"] = number_of_points;
      J["time_voronoi"] = time_voronoi;
      J["time_neighbours"] = time_neighbours;
      J["time_triangulation"] = time_triangulation;
      J["data"] = unshared_data;
      J["stack"] = unshared_stack;
      J["nb_vertices_avg"] = number_of_vertices;
      J["nb_facets_avg"] = number_of_facets;
      std::ofstream output(filename_output);
      output << std::setw(4) << J << std::endl;
      output.close();
    }
  }

  // if (number > 3) return;
  //
  // graphics::Visualizer vis;
  // vis.add_topology( transport.diagram() );
  // vis.run();

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
