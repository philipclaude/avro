#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/plots.h"

#include "mesh/facets.h"

#include "numerics/geometry.h"

#include "voronoi/optimal_transport.h"

#include "measures.h"
#include "visualize.h"

#include "json/json.hpp"

#include <fstream>
#include <iomanip>

UT_TEST_SUITE( sandbox_semidiscrete_ot_toy )

std::shared_ptr<delaunay::DensityMeasure>
get_density( const std::string& name , coord_t number )
{
  if (name == "uniform") return std::make_shared<delaunay::DensityMeasure_Uniform>(1.0);
  else if (name == "shock") return std::make_shared<DensityMeasure_Shock>(number);
  else if (name == "cone") return std::make_shared<DensityMeasure_Cone>(number);
  else if (name == "gaussian")
  {
    // gaussian
    vecd<real_t> mu(number);
    symd<real_t> sigma(number,number);
    sigma = 0;
    for (coord_t d = 0; d < number; d++)
    {
      mu(d) = 0.5;
      sigma(d,d) = 0.02;
    }
    return std::make_shared<DensityMeasure_Gaussian>(mu,sigma);
  }
  else if (name == "sphere") return std::make_shared<DensityMeasure_Sphere>(number);
  else return nullptr;
}

class DensityField : public Field<Simplex,real_t> {
public:
  DensityField( const Points& points , Topology<Simplex>& slice , const std::vector<index_t>& sites , delaunay::DensityMeasure& density ) :
    Field<Simplex,real_t>(slice,0,DISCONTINUOUS) {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<slice.nb();k++) {
      index_t s = sites[k];
      this->value(k) = density.evaluate( 0 , nullptr , points[s] );
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const {
    std::vector<std::string> result;
    result.push_back("density");
    return result;
  }
};

class MassField : public Field<Simplex,real_t> {
public:
  MassField( const Points& points , Topology<Simplex>& slice , const std::vector<index_t>& sites , const std::vector<real_t>& mass ) :
    Field<Simplex,real_t>(slice,0,DISCONTINUOUS) {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<slice.nb();k++) {
      index_t s = sites[k];
      this->value(k) = mass[s];
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const {
    std::vector<std::string> result;
    result.push_back("density");
    return result;
  }
};

UT_TEST_CASE( test1 )
{
  typedef Polytope type;

  //std::string filename = "/Users/pcaplan/Dropbox/research/publications/imr-2021-xxxx/quantization/qntz-gaussian-lbfgs-dim-4-n-10000-snapshot-iter-100.json";
  std::string filename = "/Users/pcaplan/Dropbox/research/publications/imr-2021-xxxx/quantization/qntz-cone-lbfgs-dim-4-n-10000-snapshot-iter-100.json";
  //std::string filename = "/Users/pcaplan/Dropbox/research/publications/imr-2021-xxxx/optimal_transport/sdot-uniform-dim-4-n-1000-snapshot-iter-90.json";
  //std::string filename = "/Users/pcaplan/Dropbox/research/publications/imr-2021-xxxx/performance/points-dim-4-p-6.json";


  std::fstream file;
  file.open(filename);
  avro_assert( file.is_open() );
  std::stringstream s;
  s << file.rdbuf();

  json J = json::parse( s.str().c_str() );

  std::vector<real_t> coordinates = J.at("points");

  std::vector<real_t> weights = J.at("weights");
  coord_t number = 4;//J.at("dim");
  avro_assert( number == 4 );

  coord_t dim = number+1;
  CubeDomain<type> domain(number,dim,2);

  delaunay::SemiDiscreteOptimalTransport<type> transport(domain,nullptr);
  transport.sample( coordinates.size()/number );

  transport.set_delaunay( coordinates.data() , number );
  transport.set_weights( weights.data() );

  transport.compute_laguerre();

  HyperSlice<type> slice(transport.diagram());
  delaunay::IntegrationSimplices& triangulation = transport.simplices();

  std::vector<real_t> center(number,0.01);
  slice.compute( center , 0 );

/*
  std::shared_ptr<delaunay::TriangulationCells> tc = std::make_shared<delaunay::TriangulationCells>(triangulation);
  triangulation.fields().make("c",tc);

  std::shared_ptr<delaunay::TriangulationElements> te = std::make_shared<delaunay::TriangulationElements>(triangulation);
  triangulation.fields().make("e",te);
*/

  graphics::Viewer vis;

  std::shared_ptr<SliceSites> ts = std::make_shared<SliceSites>(slice.tetrahedra(),slice.tet2site());
  slice.tetrahedra().fields().make("sites",ts);

  std::shared_ptr<delaunay::DensityMeasure> density = get_density( /*J["density"]*/ "cone" , number );
  std::shared_ptr<DensityField> dfld = std::make_shared<DensityField>( transport.delaunay() , slice.tetrahedra(),slice.tet2site() , *density.get() );
  slice.tetrahedra().fields().make("d",dfld);

  std::shared_ptr<MassField> mfld;
  if (J.contains("mass")) {
    std::vector<real_t> mass = J.at("mass");
    mfld = std::make_shared<MassField>( transport.delaunay() , slice.tetrahedra(),slice.tet2site() , mass );
    slice.tetrahedra().fields().make("m",mfld);
  }

  vis.add( slice.tetrahedra() );
  vis.add( slice.edges() );

  std::string viewfile = "";//"/Users/pcaplan/Dropbox/Codes/mach-II/avro/test/view.json";
  vis.run();
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( sandbox_semidiscrete_ot_toy )
