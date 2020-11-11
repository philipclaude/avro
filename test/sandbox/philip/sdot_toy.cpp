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

class DensityMeasure_Test : public delaunay::DensityMeasure
{
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const
  {
    //return .1 + 100*x[0]*x[0];
    return 1e1*( 1 + sin(2*M_PI*x[0])*sin(2*M_PI*x[1]) );
  }
};

class DensityMeasure_Shock : public delaunay::DensityMeasure
{
public:
  DensityMeasure_Shock( coord_t dim ) :
    dim_(dim)
  {
    k0_ = 10.0;
    k1_ = 10.0;
    vs_ = 0.7;
    r0_ = 0.4;
    alpha_ = 1.0;
  }

  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const
  {
    real_t rt = r0_ + vs_*x[dim_-1];
    real_t x2 = 0.0;
    for (coord_t d = 0; d < dim_-1; d++)
    {
      x2 += x[d]*x[d];
    }
    x2 = std::sqrt(x2);
    return 1.0 + k0_*std::exp( -alpha_*x[dim_-1] )*std::exp( -k1_*(rt-x2)*(rt-x2) );
  }
private:
  coord_t dim_;
  real_t k0_,k1_,vs_,r0_,alpha_;
};

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t number = 2;
  index_t nb_points = 1e2;

  coord_t dim = 3;
  CubeDomain<type> domain(number,dim,2);

  delaunay::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Test density2;
  DensityMeasure_Shock density3(number);

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density3);
  transport.sample( nb_points );

  transport.optimize_points(50);

  const std::vector<real_t>& mass = transport.mass();
  real_t mass_total = 0.0;
  for (index_t k = 0; k < mass.size(); k++)
    mass_total += mass[k];
  real_t mass_min = * std::min_element( mass.begin() , mass.end() );
  real_t mass_max = * std::max_element( mass.begin() , mass.end() );
  printf("total mass = %g, min = %g, max = %g, average = %g\n",mass_total,mass_min,mass_max,mass_total/real_t(nb_points));
  std::vector<real_t> nu( nb_points , mass_total / real_t(nb_points) );
  transport.set_nu( nu );
  transport.optimize_weights(300);

  HyperSlice<type> slice(transport.diagram());
  delaunay::IntegrationSimplices& triangulation = transport.simplices();
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
    //vis.add_topology(triangulation);
    vis.add_topology(point_plot);
    vis.add_topology(transport.diagram());
  }
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
