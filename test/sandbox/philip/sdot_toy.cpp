#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/optimal_transport.h"

#include "measures.h"
#include "visualize.h"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  coord_t number = 2;
  index_t nb_points = 1e2;

  coord_t dim = number+1;
  CubeDomain<type> domain(number,dim,2);

  voronoi::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Shock density2(number);
  DensityMeasure_Sphere density4(number);
  DensityMeasure_Cone density5(number);

  // gaussian
  vecd<real_t> mu(number);
  symd<real_t> sigma(number,number);
  sigma = 0;
  for (coord_t d = 0; d < number; d++)
  {
    mu(d) = 0.5;
    sigma(d,d) = 0.02;
  }
  DensityMeasure_Gaussian density3(mu,sigma);

  voronoi::SemiDiscreteOptimalTransport<type> transport(domain,&density);
  transport.save_every( 1e10 , "tmp/void" );
  transport.sample( nb_points );

  transport.weight_max() = 1e1;
  transport.quad_order() = 5;

  transport.optimize_points(100);

  const std::vector<real_t>& mass = transport.mass();

  real_t mass_total = 0.0;
  for (index_t k = 0; k < mass.size(); k++)
    mass_total += mass[k];
  real_t mass_min = * std::min_element( mass.begin() , mass.end() );
  real_t mass_max = * std::max_element( mass.begin() , mass.end() );
  printf("total mass = %g, min = %g, max = %g, average = %g\n",mass_total,mass_min,mass_max,mass_total/real_t(nb_points));
  std::vector<real_t> nu( nb_points , mass_total / real_t(nb_points) );
  transport.set_nu( nu );
  transport.optimize_weights(100);

  const std::vector<real_t>& massf = transport.mass();
  mass_total = 0.0;
  for (index_t k = 0; k < massf.size(); k++)
    mass_total += massf[k];
  mass_min = * std::min_element( massf.begin() , massf.end() );
  mass_max = * std::max_element( massf.begin() , massf.end() );
  printf("total mass = %g, min = %g, max = %g, average = %g\n",mass_total,mass_min,mass_max,mass_total/real_t(nb_points));

  HyperSlice<type> slice(transport.diagram());
  voronoi::IntegrationSimplices& triangulation = transport.simplices();
  if (number == 4)
  {
    std::vector<real_t> center(number,0.001);
    slice.compute( center , 0 );
    slice.save( "tmp/sdot-dim4-10000" );
  }
  else if (number == 3)
  {
    // export the mesh
    library::meshb writer;
    writer.open( 3 , "tmp/sdot-dim3-100000_tet.mesh" );
    writer.write( triangulation.points() );
    std::vector<index_t> refs( triangulation.nb() , 0 );
    writer.write( triangulation , refs );
    writer.close();

    // export the sites
    json J;
    J["field"] = triangulation.simplex2site();
    std::ofstream output("tmp/sdot-dim3-100000_sites.json");
    output << std::setw(4) << J << std::endl;
  }

  std::shared_ptr<voronoi::TriangulationCells> tc = std::make_shared<voronoi::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<voronoi::TriangulationElements> te = std::make_shared<voronoi::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  graphics::Viewer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());
  //for (index_t k = 0; k < point_plot.points().nb(); k++)
  //for (coord_t d = number; d < point_plot.points().dim(); d++)
  //  point_plot.points()[k][d] = 0.0;
  if (number == 4)
  {
    std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(slice.tetrahedra(),slice.tet2site());
    slice.tetrahedra().fields().make("sites",ts);
    vis.add( slice.tetrahedra() );
  }
  else
  {
    //vis.add(triangulation);
    vis.add(point_plot);
    vis.add(transport.diagram());
  }
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
