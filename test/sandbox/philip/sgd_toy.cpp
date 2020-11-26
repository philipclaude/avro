#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"

#include "visualize.h"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

class DensityMeasure_Gaussian : public delaunay::DensityMeasure
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
  coord_t number = 4;
  index_t nb_points = 2e4;

  coord_t dim = 4;
  CubeDomain<type> domain(number,dim,2);

  delaunay::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Test density2;
  DensityMeasure_Shock density3(number);

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density3);
  transport.sample( nb_points );
  transport.optimize_points(50);

  HyperSlice<type> slice(transport.diagram());

  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  if (number == 4)
  {
    std::vector<real_t> center(number,0.001);
    slice.compute( center , 0 );
    slice.save( "tmp/sdot-dim2-10000" );
  }
  else if (number == 3)
  {
    // export the mesh
    library::meshb writer;
    writer.open( 3 , "tmp/sdot-dim3-10000_tet.mesh" );
    writer.write( triangulation.points() );
    std::vector<index_t> refs( triangulation.nb() , 0 );
    writer.write( triangulation , refs );
    writer.close();

    // export the sites
    json J;
    J["field"] = triangulation.simplex2site();
    std::ofstream output("tmp/sdot-dim3-10000_sites.json");
    output << std::setw(4) << J << std::endl;
  }

  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  graphics::Visualizer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());
  if (number == 4)
  {
    std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(slice.tetrahedra(),slice.tet2site());
    slice.tetrahedra().fields().make("sites",ts);
    vis.add_topology( slice.tetrahedra() );
  }
  else
  {
    vis.add_topology(triangulation);
    //vis.add_topology(point_plot);
    //vis.add_topology(transport.diagram());
  }
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
