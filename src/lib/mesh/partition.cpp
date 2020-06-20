#include "common/process.h"
#include "common/tools.h"

#include "element/simplex.h"

#include "geometry/entity.h"

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
  local2global_ = mpi::receive<std::vector<index_t>>(sender,TAG_LOCAL2GLOBAL);

  coord_t dim = points_.dim();
  coord_t udim = points_.udim();
  index_t nb_points = x.size()/dim;
  avro_assert( u.size() == nb_points*udim );
  avro_assert( local2global_.size() == nb_points );
  avro_assert( identifiers.size() == nb_points );

  for (index_t k=0;k<nb_points;k++)
  {
    points_.create( &x[dim*k] );
    for (coord_t d=0;d<udim;d++)
      points_.u(k)[d] = u[k*udim+d];
    points_.set_entity( k , nullptr ); // TODO look up identifier
    global2local_.insert( { local2global_[k] , k } );
  }

  int min_id = *std::min_element(identifiers.begin(),identifiers.end());
  int max_id = *std::max_element(identifiers.begin(),identifiers.end());
  printf("min id = %d, max id = %d\n",min_id,max_id);
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
  index_t nb_points = local2global_.size();

  std::vector<real_t> coordinates( nb_points*dim , 0.0 );
  std::vector<real_t> parameters( nb_points*udim , 0.0 );
  std::vector<int> identifiers( nb_points );

  for (index_t k=0;k<nb_points;k++)
  {
    for (coord_t d=0;d<dim;d++)
      coordinates[k*dim+d] = points_[k][d];
    for (coord_t d=0;d<udim;d++)
      parameters[k*udim+d] = points_.u(k)[d];

    Entity* entity = points_.entity(k);
    if (entity==nullptr) identifiers[k] = -1;
    else identifiers[k] = entity->identifier();
  }

  int min_id = *std::min_element(identifiers.begin(),identifiers.end());
  int max_id = *std::max_element(identifiers.begin(),identifiers.end());
  printf("min id = %d, max id = %d\n",min_id,max_id);

  mpi::send( mpi::blocking{} , coordinates , receiver , TAG_COORDINATE );
  mpi::send( mpi::blocking{} , parameters , receiver , TAG_PARAMETER );
  mpi::send( mpi::blocking{} , identifiers , receiver , TAG_GEOMETRY );
  mpi::send( mpi::blocking{} , local2global_ , receiver , TAG_LOCAL2GLOBAL );
}

template<typename type>
void
Topology_Partition<type>::move_to_front( const std::vector<index_t>& pts )
{
  // move the points in the specified indices to the front
  // this will be useful when we need to move fixed points in an
  // adaptation chunk to the start of the points so that collapses
  // during the adaptation (which decrement indices) do not affect the
  // global maps
  avro_implement;
}

template class Partition<Simplex>;
template class Topology_Partition<Simplex>;

} // avro

#endif
