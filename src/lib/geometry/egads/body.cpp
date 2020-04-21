#include "geometry/tessellation.h"

#include "geometry/egads/body.h"
#include "geometry/egads/object.h"

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
  EGADS_ENSURE_SUCCESS( EG_getTopology( *object_ , &data_.reference , &data_.object_class ,
                        &data_.member_type , data_.data , &data_.nb_children , &data_.children , &data_.senses ) );

  printf("nb_children = %d\n",data_.nb_children);
  // loop through the children obtained from egads
  // and create the topological entities
  for (int k=0;k<data_.nb_children;k++)
  {
    // create the new entity
    std::shared_ptr<EGADS::Object> entity;
    entity = std::make_shared<EGADS::Object>(&data_.children[k],this);
    add(entity);

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

		if (pindex!=EG_indexBodyTopo(body_->object(),object_)) continue; // not this node
		index = k;
		break;

	}

	avro_assert_msg( index>=0 , "node not found" );

	topology_->add(index);
	EG_free(nodes);

}

void
Object::tessellate( BodyTessellation& tess , ego egads_tess ) const
{

  ego body_tess;
  Points& points = tess.points();

	// create the topology object
	//std::string name = objectClassName()+"-"+memberTypeName()+stringify(EG_indexBodyTopo(body_->object(),object_));

	// sorry bob, I don't like upper case
	//std::transform( name.begin() , name.end() , name.begin() , ::tolower );

	// create the topology associated with this entity
	//topology_ = smart_new(Topology<Simplex>)(points,number_);
	//topology_->setName(name);
	//topology_->setLevel( level_ );
	//topology_->setSorted( false ); // do not sort elements because EGADS provides them in a certain direction

  // get all children of the body
  std::vector<Entity*> entities;
  get_tessellatable(entities);

  for (index_t k=0;k<entities.size();k++)
  {
    std::shared_ptr<TopologyBase> topology;
    get_tessellation( entities_[k] , body_tess , tess );
  }

  for (index_t k=0;k<entities_.size();k++)
  {
    numerics::MatrixD<int> adj;
    entities_[k]->get_adjacency(adj);
  }
  /*

	// do the children
	for (index_t k=0;k<nb_children();k++)
	{
		// tessellate the child only if it isn't already done
		if (!child(k)->tessellated())
			child(k)->tessellate( tess , points );

		// now add the child's tessellation as a child to this one
		topology_->addChild( child(k)->topology() );
	}*/
}

#endif

void
Body::tessellate( BodyTessellation& tess ) const
{
  TessellationParameters& params = tess.parameters();

  // get the EGADS tessellation
  ego egads_tess;

  // retrieve the tessellation from each entity in the body
  std::map<ego,Entity_ptr>::const_iterator it;
  std::map<Entity_ptr, std::shared_ptr<TopologyBase> > topologies;
  for (it=children_.begin();it!=children_.end();++it)
  {
    ego object = it->first;
    const Entity& entity = *it->second.get();

    std::shared_ptr<TopologyBase> topology = nullptr;
    if (entity.number()==0)
    {
      //topology = get_tessellation_node( egads_tess , body_tess , params );
    }
    else if (entity.number()==1)
    {
      //get_tessellation_edge();
    }
    else if (entity.number()==2)
    {
      //get_tessellation_face();
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

    // loop through all the children of this entity
    for (index_t k=0;k<entity.nb_children();k++)
    {
      // find the topology for this child
      std::shared_ptr< TopologyBase > child_topology_base = topologies[ entity.child_smptr(k) ];

      if (params.type()=="simplex")
      {
        Topology<Simplex>& topology = static_cast<Topology<Simplex>&>(*topology_base.get());
        std::shared_ptr< Topology<Simplex> > child_topology = std::static_pointer_cast< Topology<Simplex>>(child_topology_base);
        topology.add_child(child_topology);
      }
      else if (params.type()=="quad")
      {
        printf("quads not currently supported.\n");
        avro_implement;
      }
      else
        avro_implement;
    }
  }

}


#if 0
template<>
void
Entity::retrieveTessellation<1>(ego tess)
{
  int status,idx;
  int nv,ne;
  const double *x,*u;
  index_t *edges;
  int e00,e0,e1;

  avro_assert( number_==1 );

	// degenerate edges don't have tessellations
	if (!hasTessellation())
	{
		topology_->setDummy(true);
		return;
	}

  // get the index of this object in the body
  idx = EG_indexBodyTopo( body_->object(), object_);

  // get the object tessellation from the body tessellation
  status = EG_getTessEdge( tess , idx , &nv , &x , &u );
  avro_assert( status==EGADS_SUCCESS );

  // allocate the number of edges (ne)
  ne = nv -1;
  edges = (index_t*) malloc( 2*ne*sizeof(index_t) );

  // get the index of the first point
  status = EG_localToGlobal(tess,-idx,1,&e0);
  if (status!=EGADS_SUCCESS)
    printf("could not find global index %d for edge %d (%s,%s)\n",1,-idx,
              objectClassName().c_str(),memberTypeName().c_str());
  avro_assert( status==EGADS_SUCCESS );
  e00 = e0; // needed for periodic edges

  for (int k=1;k<nv;k++)
  {
    // get the next point on the edge
    status = EG_localToGlobal(tess,-idx,k+1,&e1);
    if (status!=EGADS_SUCCESS)
       printf("could not find global index %d for edge %d (%s,%s)\n",k+1,-idx,
                objectClassName().c_str(),memberTypeName().c_str());
    avro_assert( status==EGADS_SUCCESS );

    // save the mapped indices of the edge
    edges[2*(k-1)  ] = (index_t) e0;
    edges[2*(k-1)+1] = (index_t) e1;

    // the first index is now the last one
    e0 = e1;
  }
  if(memberType_==ONENODE) edges[2*(nv-2)+1] = e00; // periodic edge

  // store the tessellation into the mesh
  std::vector<index_t> s(2);
  for (index_t k=0;k<index_t(ne);k++)
  {
		std::vector<real> uv(2,0.);
    for (index_t i=0;i<2;i++)
		{
      s[i] = edges[2*k+i] -1; // zero-bias
		}
    topology_->add(s);

  }

  free(edges);
}

template<>
void
Entity::retrieveTessellation<2>(ego tess)
{

  int status,idx;
  int nv,nt;
  const int *ptype,*pindex,*t,*ptric;
  const double *x,*u;
  index_t *triangles;
  int t0;

  avro_assert( number_==2 );

  // SHELLS, BODIES don't have explicit tessellations but their children (FACE) do
	if (!hasTessellation())
	{
		topology_->setDummy(true);
		return;
	}

  // get the index of this object in the body
  idx = EG_indexBodyTopo(body_->object(),object_);

  if (idx<0)
  {
    printf("warning: face idx = %d < 0, making positive\n",idx);
  	avro_implement;
    idx = -idx;
  }

  // get the tessellation of this object from the body tessellation
  status = EG_getTessFace( tess , idx , &nv , &x , &u , &ptype , &pindex , &nt , &t , &ptric );
  avro_assert( status==EGADS_SUCCESS );

  // allocate the triangles
  triangles = (index_t*) malloc( 3*nt*sizeof(index_t) );

  for (int k=0;k<3*nt;k++)
  {
    status = EG_localToGlobal(tess,idx,t[k],&t0);
    if (status!=EGADS_SUCCESS)
       printf("could not find global index %d for face %d (%s, %s)\n",t[k],idx,
                objectClassName().c_str(),memberTypeName().c_str());
    avro_assert( status==EGADS_SUCCESS );

    // set the triangle index to the mapped global one
    triangles[k] = (index_t) t0;
  }

  // store the tessellation into the mesh
  std::vector<index_t> s(3);
  for (index_t k=0;k<index_t(nt);k++)
  {
		std::vector<real> uv(2,0.);
    for (index_t i=0;i<3;i++)
		{
      s[i] = triangles[3*k+i] -1; // zero-bias
			//uv[0] += u[ 3*k+i ]/3.;
			//uv[1] += u[ 3*k+i ]/3.;
		}
    topology_->add(s);
  }
  free(triangles);
}

#endif

} // EGADS

} // avro
