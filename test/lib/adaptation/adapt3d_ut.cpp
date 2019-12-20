#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parameters.h"

#include "common/error.h"

#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace luna;

UT_TEST_SUITE( adaptation_adapt3d_suite)

UT_TEST_CASE(adapt_test)
{
  // setup the topology
  coord_t number = 3;

  EGADS::Context context;
  std::vector<real_t> lens(number,1.0);
  std::vector<real_t> x0(number,0.0);

  // parameters
  library::MetricField_Uniform analytic(number,0.1);

  // geometry
  EGADS::Cube box(&context,{1,1,1});

  // structured grid
  std::vector<index_t> dims(number,6);
  std::shared_ptr<Topology<Simplex>> ptopology = std::make_shared<CKF_Triangulation>(dims);
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,number);
  pmesh->add(ptopology);
  ptopology->points().copy(pmesh->points());

  pmesh->points().print();

  // tag the points onto the body
  ptopology->points().attach( box );

  // define the problem and adapt
  AdaptationParameters params;
  params.directory() = "tmp/";

  for (index_t iter=0;iter<5;iter++)
  {

    params.adapt_iter() = iter;

    // create the metric field
    std::vector<numerics::SymMatrixD<real_t>> fld;
    for (index_t k=0;k<pmesh->points().nb();k++)
      fld.push_back( analytic( pmesh->points()[k] ) );

    // create the mesh we will write to
    std::shared_ptr<Mesh> pmesh_out = std::make_shared<Mesh>(number,number);

    pmesh->retrieve<Simplex>(0).print();

    // define the problem and adapt
    AdaptationProblem problem = {*pmesh,fld,params,*pmesh_out};
    adapt<Simplex>( problem );

    // create the mesh for the next adaptation
    pmesh = std::make_shared<Mesh>(number,number);
    pmesh_out->points().copy( pmesh->points() );
    UT_ASSERT_EQUALS( pmesh->points().nb_ghost() , 0 );

    std::shared_ptr<Topology<Simplex>> ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),number);
    const Topology<Simplex>& topology_out = pmesh_out->retrieve<Simplex>(0);
    for (index_t k=0;k<pmesh_out->topology(0).nb();k++)
      ptopology->add(topology_out(k),topology_out.nv(k));
    pmesh->add(ptopology);
  }

}
UT_TEST_CASE_END(adapt_test)

UT_TEST_SUITE_END(adaptation_adapt3d_suite)
