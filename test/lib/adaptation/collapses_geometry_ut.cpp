#include "unit_tester.hpp"

#include "adaptation/adapt.h"
#include "adaptation/geometry.h"
#include "adaptation/metric.h"
#include "adaptation/parameters.h"

#include "common/error.h"

#include "geometry/egads/context.h"
#include "geometry/tessellation.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

using namespace avro;

UT_TEST_SUITE(collapses_geometry_suite)

UT_TEST_CASE(test1)
{
  // setup the topology
  coord_t number = 2;
  coord_t dim = 3;

  // parameters
  library::MetricField_Uniform analytic(2,0.1);

  // geometry
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/tire.egads");

  TessellationParameters tess_params;
  tess_params.standard();
  tess_params.min_size() = 0.2;
  tess_params.min_angle() = 20;

  ModelTessellation tess(model,tess_params);

  // create a mesh and add the topology
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(number,dim);
  pmesh->points().set_parameter_dim(number);

  // convert all the points to parameter space
  for (index_t k=0;k<tess.points().nb();k++)
  {
    std::vector<real_t> u = { tess.points().u(k,0) , tess.points().u(k,1) };
    pmesh->points().create( tess.points()[k] );
    pmesh->points().u(k,0) = u[0];
    pmesh->points().u(k,1) = u[1];
    pmesh->points().set_entity(k,tess.points().entity(k));
  }
  pmesh->points().print(true);

  // retrieve all the triangles
  Topology<Simplex> topology(pmesh->points(),2);
  tess.retrieve<Simplex>(0).get_elements( topology );
  topology.master().set_parameter(true);

  topology.close();
  topology.neighbours().compute();
  topology.inverse().build();

  Collapse<Simplex> collapser(topology);
  UT_ASSERT( collapser.master().parameter() );

  collapser.delay() = false;

  index_t pass = 0;
  while (true)
  {
    pass++;


    std::vector<bool> removed( topology.points().nb() , false );

    index_t nb_collapse = 0;

    // loop through every edge, extract the cavity and ensure every edge end point is visible
    std::vector<index_t> edges;
    topology.get_edges(edges);
    for (index_t k=0;k<edges.size()/2;k++)
    {
      index_t p0 = edges[2*k];
      index_t p1 = edges[2*k+1];

      if (removed[p0] || removed[p1]) continue;

      // skip collapses with ghost points
      if (p0<topology.points().nb_ghost() || p1<topology.points().nb_ghost())
        continue;

      // skip collapses with ghost points
      if (p0<topology.points().nb_ghost() || p1<topology.points().nb_ghost())
        continue;

      // get the entity of this edge
      Entity* entity = collapser.geometry(p0,p1);
      UT_ASSERT( entity!=nullptr );

      // skip edges along geometry Edges
      //if (entity->number()!=2) continue;

      // check if the collapse is valid
      bool accept = false;
      if (collapser.valid(p0,p1))
      {
        accept = collapser.apply( p0 , p1 );
        if (accept) removed[p0] = true;
      }
      else if (collapser.valid(p1,p0))
      {
        accept = collapser.apply( p1 , p0 );
        if (accept) removed[p1] = true;
      }
      if (accept) nb_collapse++;
    }
    if (nb_collapse==0) break;
  }

  graphics::Visualizer vis;
  vis.add_topology(topology);

  vis.run();

}
UT_TEST_CASE_END(test1)


UT_TEST_SUITE_END(collapses_geometry_suite)
