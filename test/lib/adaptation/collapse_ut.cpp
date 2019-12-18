#include "unit_tester.hpp"

#include "adaptation/primitive.h"

#include "geometry/egads/context.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/boundary.h"

using namespace luna;

UT_TEST_SUITE( collapse_test_suite )

UT_TEST_CASE(ckf_3d)
{
  coord_t number = 3;

  EGADS::Context context;

  real_t xc[3] = {0.,0.,0.};
  EGADS::Cube box(&context,{1,1,1},xc);

  std::vector<index_t> dims(number,3);
  CKF_Triangulation topology( dims );
  topology.points().print();
  topology.points().attach(box);

  bool testBoundaryVolume = false;

  topology.close();
  topology.orient();
  topology.neighbours().compute();
  topology.inverse().build();

  UT_ASSERT_EQUALS( topology.neighbours().cache().nb() , 0 );

  Collapse<Simplex> collapser(topology);

  // compute the initial volumes of the entities
  std::vector<Entity*> entities;
  //if (testBoundaryVolume)
  //  box.listTessellatableEntities(entities);
  std::vector<real_t> volume0( entities.size() );
  Boundary<Simplex> boundary0( topology );
  boundary0.extract();
  for (index_t k=0;k<entities.size();k++)
  {
    volume0[k] = boundary0.volume( boundary0.indexof( entities[k] ) );
    UT_ASSERT_NEAR( volume0[k] , 1. , 1e-12 );
  }

  index_t pass = 0;
  while (true)
  {
    pass++;

    // get the edges in the mesh
    std::vector<index_t> edges;
    topology.get_edges(edges);

    printf("performing %lu collapses\n",edges.size()/2);

    std::vector<bool> removed( topology.points().nb() , false );

    index_t nb_invalid_geometry = 0;
    index_t nb_invalid_operator = 0;

    index_t nb_collapse = 0;
    for (index_t k=0;k<edges.size()/2;k++)
    {

      index_t p0 = edges[2*k];
      index_t p1 = edges[2*k+1];

      // skip collapses with removed points
      if (removed[p0] || removed[p1]) continue;

      // skip collapses with ghost points
      if (p0<topology.points().nb_ghost() || p1<topology.points().nb_ghost())
        continue;

      // check if the collapse is valid
      bool accept = false;
      if (collapser.valid(p0,p1))
      {
        accept = collapser.apply( p0 , p1 );
        if (accept) removed[p0] = true;
        else nb_invalid_operator++;
      }
      else if (collapser.valid(p1,p0))
      {
        accept = collapser.apply( p1 , p0 );
        if (accept) removed[p1] = true;
        else nb_invalid_operator++;
      }
      else
        nb_invalid_geometry++;
      if (!accept)
        continue;

      nb_collapse++;

      if (!testBoundaryVolume)
        continue;

      // compute the boundary of the new topology
      Boundary<Simplex> boundary( topology );
      boundary.extract();

      // check the volume equals the original volume
      for (index_t j=0;j<entities.size();j++)
      {
        real_t vol = boundary.volume( boundary.indexof( entities[j] ) );
        UT_ASSERT_NEAR( vol , volume0[j] , 1e-12 );
      }

    }
    printf("collapsed %lu edges\n",nb_collapse);

    printf("nb_invalid_geometry = %lu\n",nb_invalid_geometry);
    printf("nb_invalid_operator = %lu\n",nb_invalid_operator);

    if (nb_collapse==0) break;
  }

}
UT_TEST_CASE_END(ckf_3d)

#if 0
UT_TEST_CASE(CKF_2d)
{

  exactinit(1,0,0,10,10,10);

  Context context;

  coord_t number = 2;
  typedef smart_ptr(Body) Body_ptr;

  real_t x0[3] = {5.,5.,0.};
  real_t lx = 10.,ly=10.;
  Body_ptr square = smart_new(library::EGADSSquare)(&context,x0,lx,ly);

  std::vector<real_t> lens(number,10.);
  std::vector<index_t> dims(number,30);

  library::CubeMesh mesh( lens , dims );

  mesh.points().attach( *square );


  bool testBoundaryVolume = false;

  MeshTopology<Simplex> topology( mesh.points() , mesh.number() );
  mesh.retrieveElements( number , topology );

  Cavity<Simplex> cavity(topology);

  topology.close();

  topology.orient();
  topology.neighbours().compute();
  topology.edges().compute();
  topology.inverse().build();

  UT_ASSERT_EQUALS( topology.neighbours().cache().nb() , 0 );

  Collapse<Simplex> collapser(topology);

  // compute the initial volumes of the entities
  std::vector<Entity*> entities;
  //if (testBoundaryVolume)
  //  box.listTessellatableEntities(entities);
  std::vector<real_t> volume0( entities.size() );
  Boundary<Simplex> boundary0( topology );
  boundary0.extract();
  for (index_t k=0;k<entities.size();k++)
  {
    volume0[k] = boundary0.volume( boundary0.indexof( entities[k] ) );
    UT_ASSERT_NEAR( volume0[k] , 1. , 1e-12 );
  }

  topology.points().print("v",true);

  index_t pass = 0;
  while (true)
  {
    pass++;

    // get the edges in the mesh
    std::vector<index_t> edges;
    topology.getEdges(edges);

    printf("performing %lu collapses\n",edges.size()/2);

    std::vector<bool> removed( topology.points().nb() , false );

    index_t nb_invalid_geometry = 0;
    index_t nb_invalid_operator = 0;

    index_t nb_collapse = 0;
    for (index_t k=0;k<edges.size()/2;k++)
    {

      index_t p0 = edges[2*k];
      index_t p1 = edges[2*k+1];

      // skip collapses with removed points
      if (removed[p0] || removed[p1]) continue;

      // skip collapses with ghost points
      if (p0<topology.points().nb_ghost() || p1<topology.points().nb_ghost())
        continue;

      // check if the collapse is valid
      bool accept = false;
      if (collapser.valid(p0,p1))
      {
        accept = collapser.apply( p0 , p1 );
        if (accept) removed[p0] = true;
        else nb_invalid_operator++;
      }
      else if (collapser.valid(p1,p0))
      {
        accept = collapser.apply( p1 , p0 );
        if (accept) removed[p1] = true;
        else nb_invalid_operator++;
      }
      else
        nb_invalid_geometry++;
      if (!accept)
        continue;

      nb_collapse++;

      if (!testBoundaryVolume)
        continue;

      // compute the boundary of the new topology
      Boundary<Simplex> boundary( topology );
      boundary.extract();

      // check the volume equals the original volume
      for (index_t j=0;j<entities.size();j++)
      {
        real_t vol = boundary.volume( boundary.indexof( entities[j] ) );
        UT_ASSERT_NEAR( vol , volume0[j] , 1e-12 );
      }

    }
    printf("collapsed %lu edges\n",nb_collapse);

    printf("nb_invalid_geometry = %lu\n",nb_invalid_geometry);
    printf("nb_invalid_operator = %lu\n",nb_invalid_operator);

    if (nb_collapse==0) break;
  }

  Boundary<Simplex> boundaryf(topology);
  boundaryf.extract();

  printf("nb elements = %lu\n",topology.nb());

  library::Plottable<Simplex> bplot0( boundary0 );
  library::Plottable<Simplex> bplotf( boundaryf );

  Debugger<Simplex> debugger( topology );

  GET_PLOTTER(plotter);
  plotter->addPlot(debugger);
  //plotter->addPlot(bplot0);
  plotter->addPlot(bplotf);
  plotter->run();
}
UT_TEST_CASE_END(CKF_2d)

#endif

UT_TEST_SUITE_END( collapse_test_suite )
