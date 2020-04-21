#include "geometry/tessellation.h"
#include "geometry/body.h"
#include "geometry/entity.h"
#include "geometry/model.h"

#include "numerics/geometry.h"

#include <egads.h>

namespace avro
{

ModelTessellation::ModelTessellation( Model& model , TessellationParameters& params) :
  Mesh(3,0), model_(model)
{

  // need to set the number of the mesh (it's initialized to 0)
  for (index_t k=0;k<model.nb_bodies();k++)
  {
    if (model.body(k).number()>number_)
      number_ = model.body(k).number();
  }

  // create a dummy topology to hold the children
  // which has the same topological number of the mesh
  avro_implement;
  //root_ = smart_new(Topology<Simplex>)(vertices_,number_);
  //root_->setName(model.name());
  //root_->setDummy(true);

  // set the dimension of the internal vertices to be filled
  internal_points_.set_dim( points_.dim() );

  volume_ = 0.;
  for (index_t k=0;k<model.nb_bodies();k++)
  {
    // speedup hack for wazowski demo
    //if (model.body(k)->number()!=2) continue;

    // tessellate the body in the model
    std::shared_ptr<BodyTessellation> body_tess;
    if (params.type() == "simplex")
    {
      body_tess = std::make_shared<BodyTessellation>( model.body(k) , params );
    }
    else
      avro_implement;

    /*
    int sign = ( model.is_interior(k) ) ? -1 : 1;
    real db = body_tess.volume();
    //printf("body %d, sign = %d, volume = %g\n",int(k),sign,db);
    volume_ = volume_ +sign*db;
    */

    // save the original number of vertices
    index_t nb0 = points_.nb();

    // copy the body tessellation into the model's tessellation
    copy_mesh( *body_tess.get() );

    // the vertex needs to receive information as to which body it's on
    // which is the number of children of the root topology
    for (index_t j=nb0;j<points_.nb();j++)
      points_.body(j) = k +1;

    // copy the interior/exterior flag
    //interior_.push_back( model.interior(k) );

    // retrieve the internal points from the body tessellation
    /*
    if (model.is_interior(k))
    {
      body_tess->make_internal_points();
      get_body_internal_points(bodyTess);
    }
    */

    //add(body_tess);
  }
}

void
ModelTessellation::copy_mesh( const BodyTessellation& body_tess )
{
/*
  for (index_t k=0;k<body_tess.nb_topologies();k++)
  {
    // create the new topology referencing the model tessellation's vertices
    // using the correct topological number
    //Topology_ptr t = std::make_shared<Topology<Simplex>>(points_,body_tess.topology(k).number());

    // perform a deepcopy of the topology
    //body_tess.topology(k).deepcopy( t );

    // each body references its own set of vertices
    // we need to add them to these vertices but then need to adjust
    // the topology indices by the offset which is the current number of vertices
    //t->set_offset(false);
    //t->offset_by( points_.nb() );

    // add the topology to the root (dummy) topology
    //root_->add_child(t);

  }

  // add all the vertices from the body tessellation
  for (index_t k=0;k<body_tess.points().nb();k++)
  {

    index_t id = points_.nb();

    // create the vertex coordinates
    points_.create( body_tess.points()[k] );
    points_.set_entity( id , body_tess.points().entity(k) );
  }
  */
}

Model&
ModelTessellation::model() const
{
  return model_;
}

void
ModelTessellation::get_body_internal_points( const BodyTessellation& body_tess )
{
  for (index_t k=0;k<body_tess.nb_internal();k++)
    internal_points_.create( body_tess.internal_point(k) );
  printf("retrieved internal points\n");
}

#if 0
BodyTessellation::BodyTessellation( Body& body , const real* sizes ) :
  Mesh<Simplex>(3,body.number(),body.name()), body_(body)
{
  typedef smart_ptr(Entity) Entity_ptr;

  // call egads to tessellate the body and retrieve the tessellations of the children
  ego tess;
  int npoints,status,ptype,pindex;
  double params[3];
  double box[6];
  real coordinate[3];

  // get the bounding box of the model
  EG_getBoundingBox( body.object() , box );
  double size = box[3] -box[0];
  if (size < box[4] -box[1]) size = box[4] -box[1];
  if (size < box[5] -box[2]) size = box[5] -box[2];

  // bob's magic numbers
  if (sizes)
  {
    params[0] = sizes[0];
    params[1] = sizes[1];
  }
  else
  {
    params[0] = .25*size;
    params[1] = .01*size;
  }
  params[2] = 30.;

  // ask EGADS to tessellate the body
  EGADS_ENSURE_SUCCESS( EG_makeTessBody( body.object() , params , &tess ) );

  // get the number of vertices
  EGADS_CHECK_SUCCESS( EG_statusTessBody( tess , body.pobject() , &status , &npoints ) );
  avro_assert_msg( status==1 , "egads status tessellation = %d" , status );

  // get all the children of the body
  std::vector<Entity_ptr> children;

  // loop through the children
  for (index_t k=0;k<body.nb_entities();k++)
  {
    body.entity(k)->tessellate( tess , vertices_ );
    addTopology( body.entity(k)->topology() );
  }

  // get all the body entities
  std::vector<Entity*> entities;
  body.listTessellatableEntities(entities);

  // now get the vertices from the tessellation object
  for (index_t k=0;k<index_t(npoints);k++)
  {
    EGADS_CHECK_SUCCESS( EG_getGlobal( tess , k+1 , &ptype , &pindex , coordinate ) );
    vertices_.create(coordinate);

    // lookup which entity this is
    bool found = false;
    for (index_t j=0;j<entities.size();j++)
    {
      Entity* e = entities[j];
      if (e->number()==0 && ptype==0)
      {
        index_t idx = EG_indexBodyTopo( body.object() , e->object() );
        if (int(idx)==pindex)
        {
          found = true;
          vertices_.setEntity( k , e );
          break;
        }
      }
      if (e->number()==1 && ptype>0)
      {
        index_t idx = EG_indexBodyTopo( body.object() , e->object() );
        if (int(idx)==pindex)
        {
          found = true;
          vertices_.setEntity( k , e );
          break;
        }
      }
      if (e->number()==2 && ptype<0)
      {
        index_t idx = EG_indexBodyTopo( body.object() , e->object() );
        if (int(idx)==pindex)
        {
          found = true;
          vertices_.setEntity( k , e );
          break;
        }
      }
    }
    avro_assert( found );
  }
  avro_assert( nb_topologies()==body.nb_entities() );
}

void
BodyTessellation::makeInternalPoints()
{
  // only closed bodies: solid bodies, closed shells (2-topology) or wirebodies (1-topology)
  if (!body_.closed())
  {
    printf("body is not closed: type = %s\n",body_.memberTypeName().c_str());
    return;
  }
  else printf("body is closed, making internal points\n");

  internalVertices_.setDimension( vertices_.dim() );

  // get all the n-simplices in the lower dimensional topologies
  // the assumption here is that these n-simplices bound the body
  Topology<Simplex> boundary(vertices_,number_);
  retrieveElements( number_ , boundary );

  Mesh<Simplex> volume( vertices_.dim() , number_+1 );
  if (number_==2)
  {
    TetGen tetgen(boundary,volume);
    tetgen.call("pQ");
  }
  else if (number_==1)
  {
    Triangle triangle(boundary,volume);
    triangle.call("pzYQ");
  }
  else
    avro_assert_not_reached;
  printf("called mesher\n");

  Topology<Simplex> simplices( volume.vertices() , number_+1 );
  volume.retrieveElements( number_+1 , simplices );

  // add the centroid of each volume simplex to the internal points list
  std::vector<real> centroid( vertices_.dim() );
  for (index_t k=0;k<simplices.nb();k++)
  {
    geometrics::centroid( simplices(k) , simplices.nv(k) , simplices.vertices() , centroid );
    internalVertices_.create( centroid );
  }
}

real
BodyTessellation::volume()
{
  Topology<Simplex> boundary(vertices_,number_);
  boundary.setSorted(false);
  retrieveElements(number_,boundary);
  return fabs(boundary.calculateBoundingVolume());
}
#endif


} // avro
