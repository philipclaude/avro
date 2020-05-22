#include "unit_tester.hpp"

#include "common/error.h"

#include "geometry/egads/context.h"
#include "geometry/egads/object.h"

#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/boundary.h"
#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE( mesh_parameter_suite)

UT_TEST_CASE(test1)
{
  // geometry
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");
  //EGADS::Model model(&context,"data/bunny.stp");

  TessellationParameters tess_params;
  tess_params.standard();
  ModelTessellation tess(model,tess_params);

  // create a mesh and add the topology
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(2,2);
  pmesh->points().set_parameter_dim(2);

  // convert all the points to parameter space
  std::set<Entity*> entities;
  for (index_t k=0;k<tess.points().nb();k++)
  {
    std::vector<real_t> u = { tess.points().u(k,0) , tess.points().u(k,1) };
    pmesh->points().create( u.data() );
    pmesh->points().u(k,0) = u[0];
    pmesh->points().u(k,1) = u[1];
    pmesh->points().set_entity(k,tess.points().entity(k));
    entities.insert( pmesh->points().entity(k) );
  }
  pmesh->points().print(true);

  // retrieve all the triangles
  std::shared_ptr<Topology<Simplex>> ptopology;
  ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),2);
  pmesh->add(ptopology);
  tess.retrieve<Simplex>(0).get_elements( *ptopology );
  Topology<Simplex>& topology = *ptopology.get();

  real_t vol;

  // first orient it with the nonsensical coordinates along the edges
  topology.orient();
  vol = topology.volume();
  printf("volume = %g\n",vol);
  UT_ASSERT( vol > 1e10 );

  // now tell the master that it is in parameter space, so
  // it knows that it should retrieve appropriate geometry coordinates
  topology.master().set_parameter(true);
  topology.orient();

  std::vector<real_t> volumes( topology.nb() );
  topology.get_volumes( volumes );

  for (index_t k=0;k<volumes.size();k++)
  {
    UT_ASSERT( volumes[k] >= 0 );
    if (volumes[k] < 0) printf("volume %lu = %g\n",k,volumes[k]);
  }

  // compute the volume
  vol = topology.volume();
  printf("volume = %g\n",vol);

  std::set<Entity*>::iterator it;
  real_t total_area = 0.0;
  for (it=entities.begin();it!=entities.end();it++)
  {
    EGADS::Object* obj = static_cast<EGADS::Object*>(*it);
    if (obj->number()!=2) continue;
    real_t area;
    EG_getArea( *obj->object() , NULL , &area );
    printf("area = %g\n",area);
    total_area += area;
  }
  printf("total area = %g\n",total_area);

  // compute all the edge lengths

  Boundary<Simplex> boundary(topology);
  boundary.extract();
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(mesh_parameter_suite)
