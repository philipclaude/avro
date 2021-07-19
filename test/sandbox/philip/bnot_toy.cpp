#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

class DensityMeasure_Test : public delaunay::DensityMeasure
{
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const
  {
    return 1e1*( 1 + sin(2*M_PI*x[0])*sin(2*M_PI*x[1]) );
  }
};

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t number = 2;
  index_t nb_points = 1e3;

  coord_t dim = number+1;
  CubeDomain<type> domain(number,dim,2);

  delaunay::DensityMeasure_Uniform density(1.0);

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density);
  transport.sample( nb_points );

  transport.generate_bluenoise();

  json J;
  std::vector<real_t> x( transport.delaunay().nb()*number );
  std::vector<real_t> w( transport.delaunay().nb() );
  index_t i = 0;
  for (index_t k = 0; k < transport.delaunay().nb(); k++)
  {
    for (coord_t d = 0; d < number; d++)
      x[i++] = transport.delaunay()[k][d];
  }
  J["x"] = x;
  J["w"] = transport.weights();
  J["m"] = transport.mass();

  std::ofstream output("tmp/bluenoise-dim"+std::to_string(number)+"-n"+std::to_string(transport.delaunay().nb())+".json");
  output << std::setw(4) << J << std::endl;

  if (number == 2)
  {
    std::string filename = "bluenoise-dim"+std::to_string(number)+"-n"+std::to_string(transport.delaunay().nb())+".txt";
    FILE* fid = fopen(filename.c_str(),"w");
    fprintf(fid,"%lu\n",transport.delaunay().nb());
    for (index_t k = 0; k < transport.delaunay().nb(); k++)
    {
      fprintf(fid,"%.12e %.12e\n",transport.delaunay()[k][0],transport.delaunay()[k][1]);
    }
    fclose(fid);

  }

  if (number > 3 || (nb_points >= 1e6)) return;
  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  graphics::Viewer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());

  vis.add(triangulation);
  //vis.add(point_plot);
  //vis.add(transport.diagram());
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
