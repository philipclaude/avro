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
  std::string prefix = "tmp/sdot-dim4-1000";
  library::meshb mesh( prefix+"_tet.mesh");
  std::ifstream field( prefix+"_sites.json");
  json J = json::parse(field);
  std::vector<index_t> sites = J["field"];

  Topology<Simplex>& tet = mesh.retrieve<Simplex>(0);
  std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(tet,sites);
  tet.fields().make("sites",ts);

  graphics::Visualizer vis;
  vis.add_topology(tet);
  
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
