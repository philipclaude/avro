#include "common/process.h"
#include "common/tools.h"

#include "element/simplex.h"

#include "geometry/entity.h"

#include "mesh/boundary.h"
#include "mesh/facets.h"
#include "mesh/inverse.h"
#include "mesh/partition.h"
#include "mesh/topology.h"

#include <set>

#ifdef AVRO_MPI

#include "common/mpi.hpp"

#include "mesh/mpi_tags.h"

#include <parmetis.h>

#if PARMETIS_MAJOR_VERSION == 3
#define PARM_INT idxtype
#define PARM_REAL float
#elif PARMETIS_MAJOR_VERSION == 4
#define PARM_INT idx_t
#define PARM_REAL float
#else
#error "unknown version of parmetis"
#endif

// option to use either metis or parmetis
#define USE_PARMETIS 0

namespace avro
{

template<typename type>
Partition<type>::Partition( const Topology<type>& topology ) :
  topology_(topology)
{
  // if mesh closing is needed, it should be done after partitioning
  avro_assert( !topology_.closed() );
  initialize();
}

template<typename type>
void
Partition<type>::initialize()
{
  // extract the adjacencies
  for (index_t k=0;k<topology_.nb();k++)
  {
    for (index_t j=0;j<topology_.neighbours().nfacets();j++)
    {
      int n = topology_.neighbours()(k,j);
      if (n<0) continue;      // do not create adjacency information for boundary elements
      if (n<int(k)) continue; // only create an edge once
      edges_.push_back( k );
      edges_.push_back( index_t(n) );
    }
  }

  // convert to csr
  std::vector< std::set<index_t> > node2node(topology_.nb());
  for (index_t k=0;k<edges_.size()/2;k++)
  {
    node2node[ edges_[2*k]   ].insert( edges_[2*k+1] );
    node2node[ edges_[2*k+1] ].insert( edges_[2*k] );
  }

  xadj_.resize( topology_.nb()+1 );
  xadj_[0] = 0;
  adjncy_.clear();
  std::set<index_t>::iterator it;
  for (index_t k=0;k<topology_.nb();k++)
  {
    xadj_[k+1] = xadj_[k] + node2node[k].size();
    for (it=node2node[k].begin();it!=node2node[k].end();++it)
      adjncy_.push_back(*it);
  }
  avro_assert_msg( adjncy_.size()==edges_.size() , "|adjcny| = %lu, nb_edges = %lu" , adjncy_.size() , edges_.size()/2 );
}

template<typename type>
void
Partition<type>::compute( index_t nparts0 )
{
  // setup the element-element graph
  std::vector<PARM_INT> vtxdist(nparts0+1,0);
  std::vector<PARM_INT> xadj(xadj_.begin(),xadj_.end());
  PARM_INT *pvwgt = NULL;
  std::vector<PARM_INT> adjncy(adjncy_.begin(),adjncy_.end());
  PARM_INT *padjwgt = NULL;
  PARM_INT wgtflag = 0;

  UNUSED(padjwgt);
  UNUSED(wgtflag);

  if (nparts0==1)
  {
    // if there's only one partition, then all elements
    // will receive the same coloring
    partition_.resize( topology_.nb() , 0 );
    return;
  }

  std::vector<PARM_INT> adjwgt(adjwgt_.begin(),adjwgt_.end());
  if (weighted())
  {
    padjwgt = adjwgt.data();
    wgtflag = 1;
  }
  else
  {
    padjwgt = NULL;
    wgtflag = 0;
  }

  PARM_INT bias = 0;
  PARM_INT ncon = 1;
  PARM_INT nparts = nparts0;
  PARM_INT nvtxs = topology_.nb();

  UNUSED(bias);
  UNUSED(pvwgt);
  UNUSED(nvtxs);

  std::vector<PARM_REAL> tpwgts(ncon*nparts,1./nparts);
  std::vector<PARM_REAL> ubvec(ncon,1.05);
  PARM_INT options[3];
  UNUSED(options);

  options[0] = 0;
  options[1] = 3;
  options[2] = 0;

  PARM_INT edgecut = 0;
  std::vector<PARM_INT> part(topology_.nb());

  // ask parmetis to partition the element-element graph
  //MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm comm = mpi::comm_cast( ProcessMPI::get_comm() );
  UNUSED(comm);
#if PARMETIS_MAJOR_VERSION == 4
  int result =
#endif

#if USE_PARMETIS
  ParMETIS_V3_PartKway( vtxdist.data() , xadj.data() , adjncy.data() ,
                        pvwgt, padjwgt, &wgtflag,
                        &bias , &ncon , &nparts,
                        tpwgts.data() , ubvec.data(),
                        options,
                        &edgecut,
                        part.data(),
                        &comm);
#else
  METIS_PartGraphKway( &nvtxs , &ncon , xadj.data() , adjncy.data() ,
                       NULL , NULL , NULL , &nparts , NULL , NULL , NULL , &edgecut , part.data() );
#endif

#if PARMETIS_MAJOR_VERSION == 4
  if (result!=METIS_OK) printf("error in parmetis %d\n",result);
#endif

  // save the partitioning
  partition_.resize( topology_.nb() );
  for (index_t j=0;j<part.size();j++)
    partition_[j] = part[j];

}

template<typename type>
void
Partition<type>::get( std::vector< std::shared_ptr<Topology_Partition<type>>>& parts ) const
{
  // create the partitions
  for (index_t k=0;k<parts.size();k++)
  {
    parts[k] = std::make_shared<Topology_Partition<type>>(topology_.points().dim(),topology_.points().udim(),topology_.number());
  }

  avro_assert( partition_.size() == topology_.nb() );
  for (index_t k=0;k<topology_.nb();k++)
  {
    parts[ partition_[k] ]->add( topology_(k) , topology_.nv(k) );
  }

  for (index_t k=0;k<parts.size();k++)
  {
    parts[k]->extract_points( topology_.points() );
    parts[k]->convert();

    // this ensures that all partitions only contains points they index
    avro_assert( parts[k]->all_points_accounted() );
  }
}

template<typename type>
void
Partition<type>::compute_interface( std::vector<std::set<index_t>>& elems , std::vector<std::set<index_t>>& halo ) const
{
  Facets facets( topology_ );
  facets.compute();

  // first extract all points on the interface
  std::set<index_t> pts;
  std::vector<index_t> facet( topology_.number() );
  for (index_t k=0;k<facets.nb();k++)
  {
    // if this is on the boundary of the partition, then it is not in the interface
    if (facets.boundary(k)) continue;

    // determine if the left and right partition are the same
    index_t p0 = partition_[ facets.side0(k) ];
    index_t p1 = partition_[ facets.side1(k) ];

    // if the partitions are the same, then this is not in the interface
    if (p0==p1) continue;

    facets.retrieve( k , facet );
    for (index_t j=0;j<facet.size();j++)
      pts.insert( facet[j] );
  }

  // now extract all elements with a point on the interface
  // TODO check that the neighbours are computed since these are needed for the inverse
  InverseTopology<type> inverse(topology_);
  inverse.build();
  std::set<index_t>::iterator it;
  std::vector<index_t> ball;
  for (it=pts.begin();it!=pts.end();it++)
  {
    ball.clear();
    inverse.ball( *it , ball );

    // go through every element in the ball and add it to the list of elements
    for (index_t j=0;j<ball.size();j++)
    {
      index_t elem = ball[j];
      elems[ partition_[elem] ].insert(elem);

      for (index_t i=0;i<topology_.nv(elem);i++)
      {
        if ( pts.find(topology_(elem,i)) == pts.end() ) continue; // not in the halo
        halo[ partition_[elem] ].insert( topology_(elem,i) );
      }
    }
  }
}

template<typename type>
Topology_Partition<type>::Topology_Partition( coord_t dim , coord_t udim , coord_t number ) :
  Topology<type>(points_,number),
  points_(dim,udim)
{}

template<typename type>
void
Topology_Partition<type>::extract_points( const Points& points )
{
  local2global_.clear();
  global2local_.clear();
  points_.clear();

  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<this->nv(k);j++)
    local2global_.push_back( (*this)(k,j) );
  uniquify(local2global_);

  coord_t udim = points.udim();
  for (index_t k=0;k<local2global_.size();k++)
  {
    index_t global = local2global_[k];
    points_.create( points[global] );
    for (coord_t d=0;d<udim;d++)
      points_.u(k)[d] = points.u(global)[d];
    points_.set_entity( k , points.entity(global) );
  }

  for (index_t k=0;k<local2global_.size();k++)
  {
    global2local_.insert( {local2global_.at(k),k} );
  }
}

template<typename type>
void
Topology_Partition<type>::convert()
{
  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<this->nv(k);j++)
    (*this)(k,j) = global2local_[(*this)(k,j)];
}

template<typename type>
index_t
Topology_Partition<type>::local2global( index_t k ) const
{
  avro_assert( k < points_.nb() );
  avro_assert( k < local2global_.size() );
  return local2global_[k];
}

template<typename type>
index_t
Topology_Partition<type>::global2local( index_t k ) const
{
  avro_assert( global2local_.find(k)!=global2local_.end() );
  avro_assert( local2global(global2local_.at(k)) == k );
  return global2local_.at(k);
}

template<typename type>
Entity*
Topology_Partition<type>::lookup( int identifier , int number ) const
{
  if (identifier<0 || number<0) return nullptr; // save the lookup
  for (index_t k=0;k<entities_.size();k++)
  {
    if (entities_[k]->number()==number && entities_[k]->identifier()==identifier)
      return entities_[k];
  }
  return nullptr;
}

template<typename type>
void
Topology_Partition<type>::set_entities( const std::vector<Entity*>& entities )
{
  entities_.resize( entities.size() );
  for (index_t k=0;k<entities.size();k++)
    entities_[k] = entities[k];
}

template<typename type>
void
Topology_Partition<type>::receive( mpi::communicator& comm , index_t sender )
{
  Topology<type>::clear();
  points_.clear();
  local2global_.clear();
  global2local_.clear();

  this->data_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_INDEX);
  this->first_ = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_FIRST);
  this->last_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_LAST);

  // receive the points
  std::vector<real_t> x = mpi::receive<std::vector<real_t>>(sender,TAG_COORDINATE);
  std::vector<real_t> u = mpi::receive<std::vector<real_t>>(sender,TAG_PARAMETER);
  std::vector<int> identifiers = mpi::receive<std::vector<int>>(sender,TAG_GEOMETRY);
  std::vector<int> geometry_numbers = mpi::receive<std::vector<int>>(sender,TAG_GEOMETRY_NUMBERS);
  local2global_ = mpi::receive<std::vector<index_t>>(sender,TAG_LOCAL2GLOBAL);

  coord_t dim = points_.dim();
  coord_t udim = points_.udim();
  index_t nb_points = x.size()/dim;
  avro_assert( u.size() == nb_points*udim );
  avro_assert( identifiers.size() == nb_points );

  for (index_t k=0;k<nb_points;k++)
  {
    points_.create( &x[dim*k] );
    for (coord_t d=0;d<udim;d++)
      points_.u(k)[d] = u[k*udim+d];
    Entity* entity = lookup( identifiers[k] , geometry_numbers[k] );
    points_.set_entity( k , entity );
    global2local_.insert( { local2global_[k] , k } );
  }
}

template<typename type>
void
Topology_Partition<type>::send( mpi::communicator& comm , index_t receiver ) const
{
  // send the topology
  mpi::send( mpi::blocking{} , this->data_ ,  receiver , TAG_CELL_INDEX );
  mpi::send( mpi::blocking{} , this->first_ , receiver , TAG_CELL_FIRST );
  mpi::send( mpi::blocking{} , this->last_ ,  receiver , TAG_CELL_LAST );

  // send the points
  coord_t dim = points_.dim();
  coord_t udim = points_.udim();
  index_t nb_points = points_.nb();

  std::vector<real_t> coordinates( nb_points*dim , 0.0 );
  std::vector<real_t> parameters( nb_points*udim , 0.0 );
  std::vector<int> identifiers( nb_points , -1 );
  std::vector<int> geometry_numbers( nb_points , -1 );

  for (index_t k=0;k<nb_points;k++)
  {
    for (coord_t d=0;d<dim;d++)
      coordinates[k*dim+d] = points_[k][d];
    for (coord_t d=0;d<udim;d++)
      parameters[k*udim+d] = points_.u(k)[d];

    Entity* entity = points_.entity(k);
    if (entity!=nullptr)
    {
      identifiers[k]      = entity->identifier();
      geometry_numbers[k] = entity->number();
    }
  }

  mpi::send( mpi::blocking{} , coordinates , receiver , TAG_COORDINATE );
  mpi::send( mpi::blocking{} , parameters , receiver , TAG_PARAMETER );
  mpi::send( mpi::blocking{} , identifiers , receiver , TAG_GEOMETRY );
  mpi::send( mpi::blocking{} , geometry_numbers , receiver , TAG_GEOMETRY_NUMBERS );
  mpi::send( mpi::blocking{} , local2global_ , receiver , TAG_LOCAL2GLOBAL );
}

template<typename type>
void
Topology_Partition<type>::compute_crust()
{
  avro_assert( !this->closed() );

  crust_.clear();
  halo_.clear();

  // compute the points on the boundary
  Facets facets(*this);
  facets.compute();
  std::vector<index_t> facet( this->number() );
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;

    facets.retrieve( k , facet );

    Entity* entity = BoundaryUtils::geometryFacet( this->points_ , facet.data() , facet.size() );
    if (entity!=nullptr) continue;

    for (index_t j=0;j<facet.size();j++)
      halo_.push_back( facet[j] );
  }
  uniquify(halo_);

  // compute the crust as any point in the ball of the boundary points
  std::vector<index_t> ball;
  this->neighbours().forceCompute();
  this->neighbours().compute();
  this->inverse().clear();
  this->inverse().build();
  for (index_t k=0;k<halo_.size();k++)
  {
    ball.clear();
    this->inverse().ball( {halo_[k]},ball );

    for (index_t j=0;j<ball.size();j++)
      crust_.push_back(ball[j]);
  }
  uniquify(crust_);
}

template<typename type>
void
Topology_Partition<type>::compute_mantle()
{
  // this should only be called after the adaptation
  this->inverse().clear();
  this->inverse().build();

  std::vector<index_t> ball;
  mantle_.clear();
  for (index_t k=0;k<halo_.size();k++)
  {
    ball.clear();
    this->inverse().ball( {halo_[k]},ball );

    for (index_t j=0;j<ball.size();j++)
      mantle_.push_back(ball[j]);
  }
  uniquify(mantle_);
}

template<typename type>
void
Topology_Partition<type>::compute_halo()
{
  // recompute the neighbours
  this->neighbours().forceCompute();
  this->neighbours().fromscratch()  = true;
  this->neighbours().compute();

  halo_.clear();
  for (index_t k=0;k<this->nb();k++)
  {
    for (index_t j=0;j<this->neighbours().nfacets();j++)
    {
      // retrieve the j'th neighbour of element k
      int n = this->neighbours()(k,j);

      // skip if an interior facet
      if (n>=0) continue;

      // skip if an actual geometry facet
      Entity* entity = BoundaryUtils::geometryFacet( this->points_ , (*this)(k) , this->nv(k) );
      if (entity!=nullptr) continue;

      for (index_t i=0;i<this->neighbours().nfacets();i++)
      {
        if (i==j) continue;
        halo_.push_back( (*this)(k,i) );
      }
      continue; // no need to keep checking this element's neighbours
    }
  }
  uniquify(halo_);

  // TODO: send the boundary points back to root processor
  // it needs to save the global version of these (local) indices so it knows
  // that these points need to be re-mapped to the global points before the interface is adapted

}




template class Partition<Simplex>;
template class Topology_Partition<Simplex>;

} // avro

#endif
