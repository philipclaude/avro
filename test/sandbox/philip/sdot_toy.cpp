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
    return 1e1*( 1 + sin(2*M_PI*x[0])*sin(2*M_PI*x[1]) );
  }
};

class SliceSites : public Field<Simplex,real_t>
{
public:
  SliceSites( Topology<Simplex>& slice , const std::vector<index_t>& sites ) :
    Field<Simplex,real_t>(slice,0,DISCONTINUOUS)
  {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<slice.nb();k++)
    {
      this->value(k) = sites[k];
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("sites"); return result;}
};

UT_TEST_CASE( test1 )
{
  typedef Polytope type;
  //typedef Simplex type;
  coord_t number = 4;
  index_t nb_points = 1e3;

  coord_t dim = 5;
  CubeDomain<type> domain(number,dim,2);

  delaunay::DensityMeasure_Uniform density(1.0);
  DensityMeasure_Test density2;

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,&density2);
  transport.sample( nb_points );

  transport.optimize_points(20);

  HyperSlice<type> slice(transport.diagram());

  std::vector<real_t> center(number,0.001);
  slice.compute( center , 0 );

  delaunay::IntegrationSimplices& triangulation = transport.simplices();
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);

  std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(slice.tetrahedra(),slice.tet2site());
  slice.tetrahedra().fields().make("sites",ts);

  graphics::Visualizer vis;
  library::Plot<Simplex> point_plot(transport.delaunay());

  vis.add_topology( slice.tetrahedra() );

  //vis.add_topology(triangulation);
  //vis.add_topology(point_plot);
  //vis.add_topology(transport.diagram());
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
