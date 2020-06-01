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

UT_TEST_SUITE(cavity_geometry_suite)

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
  tess_params.min_size() = 0.1;
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

  Primitive<Simplex> primitive(topology);
  UT_ASSERT( primitive.master().parameter() );

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  std::vector<index_t> edges;
  topology.get_edges(edges);
  index_t nb_cavity = 0;
  for (index_t k=0;k<edges.size()/2;k++)
  {
    index_t p = edges[2*k];
    index_t q = edges[2*k+1];

    // get the entity of this edge
    Entity* entity = primitive.geometry(p,q);
    UT_ASSERT( entity!=nullptr );

    // skip edges along geometry Edges
    if (entity->number()==1) continue;

    nb_cavity++;

    // get the cavity elements
    std::vector<index_t> C;
    topology.intersect( {p,q} , C );
    UT_ASSERT_EQUALS( C.size() , 2 );

    // get the geomery parameters of the inserted point
    // and store them in the physical coordinates
    geometry_params( entity , topology.points() , &p , 1 , topology.points()[p] );

    // check the visibility of point p in the cavity
    primitive.set_entity( entity );
    bool accept = primitive.compute( p , topology.points()[p] , C );
    UT_ASSERT(accept);

    // get the geomery parameters of the inserted point
    // and store them in the physical coordinates
    geometry_params( entity , topology.points() , &q , 1 , topology.points()[q] );

    // check the visibility of point q in the cavity
    accept = primitive.compute( q , topology.points()[q] , C );
    UT_ASSERT(accept);

    // compute the midpoint coordinates of an insertion
    std::vector<real_t> U(2);
    std::vector<real_t> X(3);

    U[0] = 0.5*( topology.points()[p][0] + topology.points()[q][0] );
    U[1] = 0.5*( topology.points()[p][1] + topology.points()[q][1] );
    entity->evaluate(U,X);
    std::vector<real_t> us = U;

    index_t idx = topology.points().nb();
    topology.points().create( X.data() );
    topology.points().set_entity( idx , entity );
    topology.points().u(idx,0) = U[0];
    topology.points().u(idx,1) = U[1];

    accept = primitive.compute( idx , U.data() , C );
    UT_ASSERT(accept);

    std::vector<real_t> params0( 2*primitive.nodes().size() , 0. );
    geometry_params( entity , topology.points() , primitive.nodes().data() , primitive.nodes().size() , params0.data() );

    Points params(2);
    real_t u0[2] = {0,0}; // dummy coordinates for ghost
    params.create( u0 );
    for (index_t k=0;k<primitive.nodes().size();k++)
    {
      params.create( &params0[2*k] );
      params.set_entity( k+1 , topology.points().entity(primitive.nodes()[k]) );

      U[0] = params0[2*k];
      U[1] = params0[2*k+1];

      // evaluate the coordinates for the orientation check
      entity->evaluate( U , X );
      for (coord_t d=0;d<3;d++)
        topology.points()[ primitive.nodes()[k] ][d] = X[d];
    }

    // check the orientation of the original cavity
    primitive.extract_geometry(entity);
    UT_ASSERT_EQUALS( primitive.geometry().nb() , 6 ); // two real triangles + 4 ghosts (one for every edge)

    GeometryOrientationChecker checker( topology.points() , params , primitive.u2v() , entity  );
    int s = checker.signof( primitive.geometry() );
    UT_ASSERT( s>0 );

    for (index_t k=0;k<primitive.geometry().points().nb();k++)
    {
      primitive.geometry().points().set_entity( k , topology.points().entity( primitive.u2v()[k] ) );
    }
    idx = primitive.geometry().points().nb();
    primitive.geometry().points().create( us.data() );
    primitive.geometry().points().set_entity( idx , entity );

    // check the inserted point
    primitive.gcavity().set_entity( entity );
    std::vector<index_t> S( primitive.geometry().nb());
    for (index_t k=0;k<primitive.geometry().nb();k++)
      S[k] = k;

    primitive.gcavity().sign() = entity->sign(); // not actually necessary because master.parameter() will trigger the sign to be used in get_volume
    accept = primitive.gcavity().compute( idx , us.data() , S );
    UT_ASSERT(accept);

    // delete the created point (otherwise the next test will fail)
    topology.points().remove( topology.points().nb()-1 );
  }
  UT_ASSERT(nb_cavity>0);
  index_t nb_cavity0 = nb_cavity;

  // loop through every point, extract the cavity and ensure the point is visible
  for (index_t k=0;k<topology.points().nb();k++)
  {

    if (k<topology.points().nb_ghost()) continue;
    nb_cavity++;

    // get the entity of this edge
    Entity* entity = topology.points().entity(k);
    UT_ASSERT( entity!=nullptr );

    // skip edges along geometry Edges
    if (entity->number()<2) continue;

    // get the cavity elements
    std::vector<index_t> C;
    topology.intersect( {k} , C );
    index_t n0 = C.size();

    // get the geometry parameters of the inserted point
    // and store them in the physical coordinates
    geometry_params( entity , topology.points() , &k , 1 , topology.points()[k] );

    // check the visibility of point p in the cavity
    primitive.set_entity( entity );
    bool accept = primitive.compute( k , topology.points()[k] , C );
    UT_ASSERT(accept);

    // check the orientation of the original cavity
    primitive.extract_geometry(entity);
    UT_ASSERT_EQUALS( primitive.geometry().nb() , 2*n0 ); // one ghost added for every original triangle

    // compute the parameter coordinates
    std::vector<real_t> params0( 2*primitive.nodes().size() , 0. );
    geometry_params( entity , topology.points() , primitive.nodes().data() , primitive.nodes().size() , params0.data() );

    Points params(2);
    real_t u0[2] = {0,0}; // dummy coordinates for ghost
    params.create( u0 );
    std::vector<real_t> U(2),X(3);
    for (index_t k=0;k<primitive.nodes().size();k++)
    {
      params.create( &params0[2*k] );
      params.set_entity( k+1 , topology.points().entity(primitive.nodes()[k]) );

      U[0] = params0[2*k];
      U[1] = params0[2*k+1];

      // evaluate the coordinates for the orientation check
      entity->evaluate( U , X );
      for (coord_t d=0;d<3;d++)
        topology.points()[ primitive.nodes()[k] ][d] = X[d];
    }

    GeometryOrientationChecker checker( topology.points() , params , primitive.u2v() , entity  );
    int s = checker.signof( primitive.geometry() );
    UT_ASSERT( s>0 );

    for (index_t k=0;k<primitive.geometry().points().nb();k++)
    {
      primitive.geometry().points().set_entity( k , topology.points().entity( primitive.u2v()[k] ) );
    }
    index_t idx = primitive.geometry().points().nb() -1;

    primitive.gcavity().set_entity( entity );
    std::vector<index_t> S( primitive.geometry().nb());
    for (index_t k=0;k<primitive.geometry().nb();k++)
      S[k] = k;

    primitive.gcavity().sign() = entity->sign(); // not actually necessary because master.parameter() will trigger the sign to be used in get_volume
    accept = primitive.gcavity().compute( idx , topology.points().u(k) , S );
    UT_ASSERT(accept);
  }
  printf("total cavities considered = %lu\n",nb_cavity);
  UT_ASSERT(nb_cavity>nb_cavity0);
}
UT_TEST_CASE_END(test1)


UT_TEST_CASE(test2)
{
  // setup the topology
  coord_t number = 2;
  coord_t dim = 3;

  // parameters
  library::MetricField_Uniform analytic(2,0.1);

  // geometry
  EGADS::Context context;
  EGADS::Model model(&context,BASE_TEST_DIR+"/geometry/cube-cylinder.egads");

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

  Primitive<Simplex> primitive(topology);
  UT_ASSERT( primitive.master().parameter() );

  SurfaceCavity<Simplex> surface(topology);
  UT_ASSERT( surface.master().parameter() );

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  std::vector<index_t> edges;
  topology.get_edges(edges);
  index_t nb_cavity = 0;
  for (index_t k=0;k<edges.size()/2;k++)
  {
    index_t p = edges[2*k];
    index_t q = edges[2*k+1];

    // get the entity of this edge
    Entity* entity = primitive.geometry(p,q);
    UT_ASSERT( entity!=nullptr );

    // skip edges along geometry Edges
    if (entity->number()==1) continue;

    nb_cavity++;

    // get the cavity elements
    std::vector<index_t> C;
    topology.intersect( {p,q} , C );
    UT_ASSERT_EQUALS( C.size() , 2 );

    surface.extract( C , entity );

    bool accept;

    accept = surface.visible( p );
    UT_ASSERT(accept);

    accept = surface.visible( q );
    UT_ASSERT(accept);

    // compute the midpoint coordinates of an insertion
    std::vector<real_t> U(2);
    std::vector<real_t> X(3);

    // revert the coordinates back to parameter space
    geometry_params( entity , topology.points() , &p , 1 , topology.points()[p] );
    geometry_params( entity , topology.points() , &q , 1 , topology.points()[q] );

    U[0] = 0.5*( topology.points()[p][0] + topology.points()[q][0] );
    U[1] = 0.5*( topology.points()[p][1] + topology.points()[q][1] );
    entity->evaluate(U,X);

    index_t idx = topology.points().nb();
    topology.points().create( X.data() );
    topology.points().set_entity( idx , entity );
    topology.points().u(idx,0) = U[0];
    topology.points().u(idx,1) = U[1];

    accept = surface.visible( idx );
    UT_ASSERT(accept);

    // check the normals of the original cavity
    accept = surface.check_normals();
    UT_ASSERT(accept);

    // check the normals of the updated topology
    accept = surface.cavity_visible( idx );
    UT_ASSERT( accept );

    accept = surface.check_normals();
    UT_ASSERT(accept);

    // delete the created point (otherwise the next test will fail)
    topology.points().remove( topology.points().nb()-1 );
  }
  UT_ASSERT(nb_cavity>0);
  index_t nb_cavity0 = nb_cavity;

  // loop through every edge, extract the cavity and ensure every edge end point is visible
  for (index_t k=0;k<topology.points().nb();k++)
  {
    if (k<topology.points().nb_ghost()) continue;

    // get the entity of this edge
    Entity* entity = topology.points().entity(k);
    UT_ASSERT( entity!=nullptr );

    // skip edges along geometry Edges
    if (entity->number()<2) continue;

    nb_cavity++;

    // get the cavity elements
    std::vector<index_t> C;
    topology.intersect( {k} , C );

    surface.extract( C , entity );

    bool accept;

    accept = surface.visible( k );
    UT_ASSERT(accept);

    // check the normals of the original cavity
    accept = surface.check_normals();
    UT_ASSERT(accept);
  }
  UT_ASSERT(nb_cavity>nb_cavity0);

  printf("total cavities considered = %lu\n",nb_cavity);

}
UT_TEST_CASE_END(test2)


UT_TEST_SUITE_END(cavity_geometry_suite)
