//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/tessellation.h"

#include "geometry/egads/body.h"
#include "geometry/egads/data.h"
#include "geometry/egads/object.h"

#include "numerics/geometry.h"

#include <egads.h>

namespace avro
{

namespace EGADS
{

void
Body::build_hierarchy()
{
  avro_assert( nb_entities()==0 );

  // get the topology of the associated ego with the children
  EGADS_ENSURE_SUCCESS( EG_getTopology( object_ , &data_.reference , &data_.object_class ,
                        &data_.member_type , data_.data , &data_.nb_children , &data_.children , &data_.senses ) );

  // loop through the children obtained from egads
  // and create the topological entities
  for (int k=0;k<data_.nb_children;k++)
  {
    // create the new entity
    std::shared_ptr<EGADS::Object> entity;
    entity = std::make_shared<EGADS::Object>(&data_.children[k],this);
    add(entity);
    add_child( data_.children[k] , entity );

    // build the entity hierarchy
    entity->build_hierarchy();
    entity->set_parent(NULL);
    entity->set_body(this);
  }

  // determine the topological number from the children
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);

  // assign the parent hierarchy
  build_parents();
}

Body::Entity_ptr
Body::child(index_t k)
{
  avro_assert( k < nb_entities() );
  return entity_[k];
}

void
Body::add_child( ego object , Entity_ptr entity )
{
  children_.insert( {object,entity} );
}

Body::Entity_ptr
Body::lookup( ego object ) const
{
  if (children_.find(object)==children_.end())
    return nullptr;
  return children_.at(object);
}

void
Body::print() const
{
  printf("EGADS: number = %u , class = %s, type = %s at %p\n",number_,
  EGADS::utilities::object_class_name(data_.object_class).c_str(),
  EGADS::utilities::member_type_name(data_.object_class,data_.member_type).c_str(),
  (void*)(this) );

  for (index_t k=0;k<nb_entities();k++)
  {
    printf("\t");
    entity_[k]->print();
  }
}

#if 0
void
get_tessellation_node( Entity* node , ego tess , BodyTessellation& body_tess )
{
	int status,npoints,ptype,pindex;
	real_t coordinate[3];
	ego* nodes;
	int nbNodes;
	int index = -1;

	// get the number of vertices
	EGADS_CHECK_SUCCESS( EG_statusTessBody( tess , body_->pobject() , &status , &npoints ) );
	avro_assert_msg( status==1 , "egads status tessellation = %d" , status );

	// get all nodes
  EGADS_CHECK_SUCCESS( EG_getBodyTopos( body_tess.body().object() , NULL , NODE , &nbNodes , &nodes) );

	// look for the matching ego of this node
	for (index_t k=0;k<index_t(npoints);k++)
	{
		// get global information about this vertex
		EGADS_CHECK_SUCCESS( EG_getGlobal( tess , k+1 , &ptype , &pindex , coordinate ) );

		if (ptype!=0) continue; // not a node vertex

		if (pindex!=EG_indexBodyTopo(body_->object(),&object_)) continue; // not this node
		index = k;
		break;

	}

	avro_assert_msg( index>=0 , "node not found" );

	topology_->add(index);
	EG_free(nodes);

}

#endif

std::shared_ptr<TopologyBase>
get_tessellation_node( ego body , ego egads_tess , ego node , BodyTessellation& tess , TessellationParameters& params , index_t offset )
{
  // allocate the appropriate object to hold the tessellation
  std::shared_ptr<TopologyBase> topology = nullptr;
  std::string type = params["type"];
  if (type == "simplex")
    topology = std::make_shared< Topology<Simplex> >( tess.model_points() , 0 );
  else
    avro_implement;

  int status,nb_points,ptype,pindex;
	real_t coordinate[3];
	index_t index;

  // get the number of points in the tessellation
  EGADS_CHECK_SUCCESS( EG_statusTessBody( egads_tess , &body , &status , &nb_points ) );
  avro_assert_msg( status==1 , "egads status tessellation = %d" , status );

	// look for the matching ego of this node
  index = nb_points;
	for (index_t k=0;k<index_t(nb_points);k++)
	{
		// get global information about this vertex
		EGADS_CHECK_SUCCESS( EG_getGlobal( egads_tess , k+1 , &ptype , &pindex , coordinate ) );

		if (ptype!=0) continue; // not a node

		if (pindex!=EG_indexBodyTopo(body,node)) continue; // not this node
		index = k;
		break;
	}
  avro_assert( index < index_t(nb_points) );

  // add the element
  index += offset;
  topology->add( &index , 1 );

  return topology;
}

std::shared_ptr<TopologyBase>
get_tessellation_edge( ego body , ego egads_tess , const Object& edge , BodyTessellation& tess , TessellationParameters& params , index_t offset )
{
  // allocate the appropriate object to hold the tessellation
  std::shared_ptr<TopologyBase> topology = nullptr;
  std::string type = params["type"];
  if (type == "simplex")
    topology = std::make_shared< Topology<Simplex> >( tess.model_points() , 1 );
  else
    avro_implement;

  if (!edge.tessellatable()) return topology;

  int status,idx;
  int nv,ne;
  const double *x,*u;
  index_t *edges;
  int e00,e0,e1;

  // get the index of this object in the body
  idx = EG_indexBodyTopo( body , *edge.object() );

  // get the object tessellation from the body tessellation
  status = EG_getTessEdge( egads_tess , idx , &nv , &x , &u );
  avro_assert( status==EGADS_SUCCESS );

  // allocate the number of edges (ne)
  ne = nv -1;
  edges = (index_t*) malloc( 2*ne*sizeof(index_t) );

  // get the index of the first point
  status = EG_localToGlobal(egads_tess,-idx,1,&e0);
  if (status!=EGADS_SUCCESS)
    printf("could not find global index %d for edge %d: member_type = %s, object_class = %s\n",
      1,-idx,utilities::member_type_name(edge.object_class(),edge.member_type()).c_str(),utilities::object_class_name(edge.object_class()).c_str());
  avro_assert( status==EGADS_SUCCESS );
  e00 = e0; // needed for periodic edges

  tess.model_points().set_param( e0-1 , {u[0],1e20} );
  for (int k=1;k<nv;k++)
  {
    // get the next point on the edge
    status = EG_localToGlobal(egads_tess,-idx,k+1,&e1);
    if (status!=EGADS_SUCCESS)
       printf("could not find global index %d for edge %d\n",k+1,-idx);
    avro_assert( status==EGADS_SUCCESS );

    // save the mapped indices of the edge
    edges[2*(k-1)  ] = (index_t) e0 -1;
    edges[2*(k-1)+1] = (index_t) e1 -1;

    tess.model_points().set_param( e1-1 + offset , {u[k],1e20} );

    // the first index is now the last one
    e0 = e1;
  }
  if(edge.member_type()==ONENODE) edges[2*(nv-2)+1] = e00; // periodic edge

  // store the tessellation into the mesh
  std::vector<index_t> s(2);
  for (index_t k=0;k<index_t(ne);k++)
  {
    for (index_t i=0;i<2;i++)
		{
      s[i] = edges[2*k+i] + offset;
		}
    topology->add(s.data(),s.size());
  }

  free(edges);

  return topology;
}

std::shared_ptr<TopologyBase>
get_tessellation_face( ego body , ego egads_tess , const Object& face , BodyTessellation& tess , TessellationParameters& params , index_t offset )
{
  // allocate the appropriate object to hold the tessellation
  std::shared_ptr<TopologyBase> topology = nullptr;
  std::string type = params["type"];
  if (type == "simplex")
    topology = std::make_shared< Topology<Simplex> >( tess.model_points() , 2 );
  else
    avro_implement;

  if (!face.tessellatable()) return topology;

  int status,idx;
  int nv,nt;
  const int *ptype,*pindex,*t,*ptric;
  const double *x,*u;
  index_t *triangles;
  int t0;

  // get the index of this object in the body
  idx = EG_indexBodyTopo(body,*face.object());

  if (idx<0)
  {
    printf("warning: face idx = %d < 0, making positive\n",idx);
  	avro_implement;
    idx = -idx;
  }

  // get the tessellation of this object from the body tessellation
  status = EG_getTessFace( egads_tess , idx , &nv , &x , &u , &ptype , &pindex , &nt , &t , &ptric );
  avro_assert( status==EGADS_SUCCESS );

  // allocate the triangles
  triangles = (index_t*) malloc( 3*nt*sizeof(index_t) );

  for (int k=0;k<3*nt;k++)
  {
    status = EG_localToGlobal(egads_tess,idx,t[k],&t0);
    if (status!=EGADS_SUCCESS)
       printf("could not find global index %d for face %d\n",t[k],idx);
    avro_assert( status==EGADS_SUCCESS );

    // set the triangle index to the mapped global one
    triangles[k] = (index_t) t0 -1; // zero-bias
  }

  // store the tessellation into the mesh
  std::vector<index_t> s(3);
  for (index_t k=0;k<index_t(nt);k++)
  {
    for (index_t i=0;i<3;i++)
		{
      s[i] = triangles[3*k+i] + offset;

      // save the parameter coordinates
      index_t v = t[3*k+i] -1;
      if (ptype[v]!=-1) continue; // not internal
      tess.model_points().set_param( s[i] , &u[2*v] );
		}
    topology->add(s.data(),s.size());
  }

  free(triangles);

  return topology;
}

void
Body::tessellate( BodyTessellation& tess ) const
{
  TessellationParameters& params = tess.parameters();

  // call egads to tessellate the body and retrieve the tessellations of the children
  ego egads_tess;
  int nb_points,status;
  double sizes[3];
  double box[6],size;

  // get the bounding box of the model
  EG_getBoundingBox( object_ , box );
  size = box[3] -box[0];
  if (size < box[4] -box[1]) size = box[4] -box[1];
  if (size < box[5] -box[2]) size = box[5] -box[2];

  real_t min_size = params["min size"];
  real_t min_length = params["min length"];

  sizes[0] = min_size*size;
  sizes[1] = min_length*size;
  sizes[2] = params["min angle"];

  printf("sizes = (%g,%g,%g)\n",sizes[0],sizes[1],sizes[2]);

  // ask EGADS to tessellate the body
  EGADS_ENSURE_SUCCESS( EG_makeTessBody( object_ , sizes , &egads_tess ) );

  // get the number of points
  EGADS_CHECK_SUCCESS( EG_statusTessBody( egads_tess , const_cast<ego*>(&object_) , &status , &nb_points ) );
  avro_assert_msg( status==1 , "egads status tessellation = %d" , status );

  std::shared_ptr<TopologyBase> root;
  std::string type = params["type"];
  if (type == "simplex")
    root = std::make_shared<Topology<Simplex>>(tess.model_points(),tess.number());
  else
    avro_implement;
  tess.add(root);

  index_t offset = tess.model_points().nb();

  // add all points to the tessellation
  std::map<ego,Entity_ptr>::const_iterator it;
  for (index_t k=0;k<index_t(nb_points);k++)
  {
    int ptype,pindex;
    double coordinate[3];
    EGADS_CHECK_SUCCESS( EG_getGlobal( egads_tess , k+1 , &ptype , &pindex , coordinate ) );
    tess.model_points().create(coordinate);

    // lookup which entity this is
    bool found = false;
    for (it=children_.begin();it!=children_.end();++it)
    {
      Entity* e = it->second.get();
      ego egads_object = it->first;

      if (!e->tessellatable()) continue; // don't assign non-tessellatable entities (e.g. loops)

      if (e->number()==0 && ptype==0)
      {
        index_t idx = EG_indexBodyTopo( object_ , egads_object );
        if (int(idx)==pindex)
        {
          found = true;
          tess.model_points().set_entity( k , e );
          break;
        }
      }
      if (e->number()==1 && ptype>0)
      {
        index_t idx = EG_indexBodyTopo( object_ , egads_object );
        if (int(idx)==pindex)
        {
          found = true;
          tess.model_points().set_entity( k , e );
          break;
        }
      }
      if (e->number()==2 && ptype<0)
      {
        index_t idx = EG_indexBodyTopo( object_ , egads_object );
        if (int(idx)==pindex)
        {
          found = true;
          tess.model_points().set_entity( k , e );
          break;
        }
      }
    }
    if (!found)
    {
      print();
      for (it=children_.begin();it!=children_.end();++it)
      {
        if (!it->second->tessellatable()) continue;
        it->second->print();
      }
    }
    avro_assert_msg( found , "could not find entity for pindex = %d, ptype = %d" , pindex , ptype );

  } // loop over points

  // retrieve the tessellation from each entity in the body
  std::map<Entity_ptr, std::shared_ptr<TopologyBase> > topologies;
  for (it=children_.begin();it!=children_.end();++it)
  {
    const Entity& entity = *it->second.get();
    ego egads_entity = it->first;

    std::shared_ptr<TopologyBase> topology = nullptr;
    if (entity.number()==0)
    {
      topology = get_tessellation_node( object_ , egads_tess , egads_entity , tess , params , offset );
    }
    else if (entity.number()==1)
    {
      const Object& edge = static_cast<const Object&>(entity);
      topology = get_tessellation_edge( object_ , egads_tess , edge , tess , params , offset );
    }
    else if (entity.number()==2)
    {
      const Object& face = static_cast<const Object&>(entity);
      topology = get_tessellation_face( object_ , egads_tess , face , tess , params , offset );
    }
    else
      avro_assert_not_reached;

    topologies.insert( {it->second , topology} );
  }

  // loop through the entities and construct the children of the topologies
  std::map<Entity_ptr, std::shared_ptr<TopologyBase> >::iterator itt;
  for (itt=topologies.begin();itt!=topologies.end();++itt)
  {
    const Entity& entity = *itt->first.get();
    std::shared_ptr< TopologyBase > topology_base = itt->second;
    topology_base->offset_by( tess.points().nb() );

    // if the child is of the appropriate topological number, add it to the root topology
    if (topology_base->number()==tess.number())
      root->add_child_type(topology_base);

    // loop through all the children of this entity
    for (index_t k=0;k<entity.nb_children();k++)
    {
      // find the topology for this child
      std::shared_ptr< TopologyBase > child_topology_base = topologies[ entity.child_smptr(k) ];

      if (type == "simplex")
      {
        Topology<Simplex>& topology = static_cast<Topology<Simplex>&>(*topology_base.get());
        std::shared_ptr< Topology<Simplex> > child_topology = std::static_pointer_cast< Topology<Simplex>>(child_topology_base);
        topology.add_child(child_topology);
      }
      else if (type == "quad")
      {
        printf("quads not currently supported.\n");
        avro_implement;
      }
      else
        avro_implement;
    }
  }

  for (index_t k=0;k<tess.points().nb();k++)
  {
    std::vector<real_t> u(tess.points().u(k) , tess.points().u(k)+2 );
    std::vector<real_t> x(3);
    tess.points().entity(k)->evaluate( u , x );

    real_t d = numerics::distance( x.data() , tess.points()[k] , 3 );

    // todo get actual tolerance from egads object...wait until this check fails until implementing this
    avro_assert( d < 1e-12 );

  }
}

} // EGADS

} // avro
