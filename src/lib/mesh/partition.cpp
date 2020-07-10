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
  avro_assert_msg( adjncy_.size()==edges_.size() , "|adjncy| = %lu, nb_edges = %lu" , adjncy_.size() , edges_.size()/2 );
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
  PARM_INT options[METIS_NOPTIONS];
  UNUSED(options);

  #if USE_PARMETIS
  options[0] = 0;
  options[1] = 3;
  options[2] = 0; // default random seed
  #else
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_SEED] = 0;
  #endif

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
                       NULL , NULL , NULL , &nparts , NULL , NULL , options , &edgecut , part.data() );
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
    parts[k]->convert( topology_.points() );
    avro_assert( parts[k]->all_points_accounted() );
  }
}

template<typename type>
void
Partition<type>::compute_interface_points( std::vector<std::set<index_t>>& halo ) const
{
  Facets facets( topology_ );
  facets.compute();

  // first extract all points on the interface
  std::vector<index_t> facet( topology_.number() );
  for (index_t k=0;k<facets.nb();k++)
  {
    // if this is on the boundary of the topology, then it is not in the interface
    if (facets.boundary(k)) continue;

    // determine if the left and right partition are the same
    index_t p0 = partition_[ facets.side0(k) ];
    index_t p1 = partition_[ facets.side1(k) ];

    // if the partitions are the same, then this is not in the interface
    if (p0==p1) continue;

    // retrieve the indices of the facet
    facets.retrieve( k , facet );

    if (fixed_facet(facet,topology_.points()))
    {
      // if all points in the facet are fixed, then this is not a partition boundary
      continue;
    }

    // associate any point in the facet with both partitions
    for (index_t j=0;j<facet.size();j++)
    {
      halo[p0].insert( facet[j] );
      halo[p1].insert( facet[j] );
    }
  }
}

template<typename type>
void
Partition<type>::print() const
{
  for (index_t k=0;k<topology_.nb();k++)
  {
    print_inline( topology_.get(k) , "element " + std::to_string(k) + " on partition " + std::to_string(partition_[k]) );
  }
}

template<typename type>
Topology_Partition<type>::Topology_Partition( coord_t dim , coord_t udim , coord_t number ) :
  Topology<type>(points_,number),
  points_(dim,udim)
{}

template<typename type>
void
Topology_Partition<type>::convert( const Points& points )
{
  points_.clear();

  std::vector<index_t> local2global;
  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<this->nv(k);j++)
    local2global.push_back( (*this)(k,j) );
  uniquify(local2global);

  std::map<index_t,index_t> global2local;
  for (index_t k=0;k<local2global.size();k++)
  {
    index_t global = local2global[k];
    global2local.insert({global,k});
    points_.create( points[global] );
    points_.set_param( k , points.u(global) );
    points_.set_entity( k , points.entity(global) );
    points_.set_fixed( k , points.fixed(global) );
    points_.set_global( k , points.global(global) );
  }

  for (index_t k=0;k<this->nb();k++)
  for (index_t j=0;j<this->nv(k);j++)
    (*this)(k,j) = global2local[(*this)(k,j)];
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
Topology_Partition<type>::receive( index_t sender )
{
  Topology<type>::clear();
  points_.clear();

  this->data_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_INDEX);
  this->first_ = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_FIRST);
  this->last_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_LAST);

  // receive the points
  std::vector<real_t> x = mpi::receive<std::vector<real_t>>(sender,TAG_COORDINATE);
  std::vector<real_t> u = mpi::receive<std::vector<real_t>>(sender,TAG_PARAMETER);
  std::vector<int> identifiers = mpi::receive<std::vector<int>>(sender,TAG_GEOMETRY);
  std::vector<int> geometry_numbers = mpi::receive<std::vector<int>>(sender,TAG_GEOMETRY_NUMBERS);
  std::vector<int> global = mpi::receive<std::vector<int>>(sender,TAG_LOCAL2GLOBAL);
  std::vector<int> fixed = mpi::receive<std::vector<int>>(sender,TAG_FIXED);

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
    if (fixed[k]==1)
      points_.set_fixed(k,true);
    else
      points_.set_fixed(k,false);
    points_.set_global( k , global[k] );
  }
}

template<typename type>
void
Topology_Partition<type>::send( index_t receiver ) const
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
  std::vector<int> fixed( nb_points , 0 );
  std::vector<int> global( nb_points , -1 );

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
    if (points_.fixed(k)) fixed[k] = 1;
    global[k] = points_.global(k);
  }

  mpi::send( mpi::blocking{} , coordinates , receiver , TAG_COORDINATE );
  mpi::send( mpi::blocking{} , parameters , receiver , TAG_PARAMETER );
  mpi::send( mpi::blocking{} , identifiers , receiver , TAG_GEOMETRY );
  mpi::send( mpi::blocking{} , geometry_numbers , receiver , TAG_GEOMETRY_NUMBERS );
  mpi::send( mpi::blocking{} , global , receiver , TAG_LOCAL2GLOBAL );
  mpi::send( mpi::blocking{} , fixed , receiver , TAG_FIXED );
}

template<typename type>
bool
Topology_Partition<type>::check( const Points& points ) const
{
  bool result = true;
  for (index_t j=0;j<points_.nb();j++)
  {
    real_t d = numerics::distance( points_[j] , points[points.global(j)] , points_.dim() );
    if (d > 1e-12)
    {
      printf("local point %lu -> global point %lu do not match d = %g.\n",j,points.global(j),d);
      points_.print(j,true);
      points.print(points.global(j),true);
      result = false;
    }
  }
  return result;
}

template<typename type>
PartitionBoundary<type>::PartitionBoundary( coord_t number ) :
  Topology<type>(dummy_,number),
  dummy_(number)
{
  this->set_sorted(true);
}

template<typename type>
void
PartitionBoundary<type>::compute( const Topology<type>& topology , index_t partition )
{
  Facets facets(topology);
  facets.compute();

  std::vector<index_t> facet(topology.number());
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,facet);
    Entity* entity = BoundaryUtils::geometryFacet( topology.points() , facet.data() , facet.size() );
    if (entity!=nullptr) continue;

    index_t elemL = facets.side0(k);
    index_t elemL_global = elemL; // TODO keep tracking of global element number too

    // add the global indices
    for (index_t i=0;i<facet.size();i++)
      facet[i] = topology.points().global(facet[i]);

    std::sort( facet.begin() , facet.end() );
    ElementIndices f;
    f.dim = this->number();
    f.indices = facet;

    add_left( f , elemL , elemL_global , partition );
  }
}

template<typename type>
void
PartitionBoundary<type>::add_left( const ElementIndices& f , index_t elemL , index_t elemL_global , index_t partL )
{
  ElementIndices facet;
  facet.dim = f.dim;
  facet.indices = f.indices;
  facet_.insert( {facet,facet_.size()} );

  elemL_.push_back(elemL);
  elemL_global_.push_back(elemL_global);
  partL_.push_back(partL);

  elemR_.push_back(0);
  elemR_global_.push_back(0);
  partR_.push_back(0);
}

template<typename type>
void
PartitionBoundary<type>::append( const PartitionBoundary<type>& boundary )
{
  index_t elemL,partL;

  const std::map<ElementIndices,index_t>& facets = boundary.facets();
  for (std::map<ElementIndices,index_t>::const_iterator it=facets.begin();it!=facets.end();++it)
  {
    // make sure this facet does not already exist
    const ElementIndices& facet = it->first;
    index_t k = it->second;
    int id = find(facet);
    if (id>=0)
    {
      // add information to the right
      elemR_[id] = boundary.elemL(k);
      partR_[id] = boundary.partL(k);
      continue;
    }

    elemL = boundary.elemL(k);
    partL = boundary.partL(k);

    // add the facet to this partition boundary
    add_left( facet , elemL , elemL , partL );
  }
}

template<typename type>
void
PartitionBoundary<type>::fill( const PartitionBoundary<type>& interface )
{
  for (std::map<ElementIndices,index_t>::iterator it=facet_.begin();it!=facet_.end();++it)
  {
    const ElementIndices& facet = it->first;
    index_t k = it->second;

    int id = interface.find(facet);
    avro_assert(id>=0);

    // check if we match the left element
    index_t elemL = interface.elemL(id);
    index_t partL = interface.partL(id);
    if (elemL_[k] == elemL && partL_[k] == partL)
    {
      elemR_[k] = interface.elemR(id);
      partR_[k] = interface.partR(id);
    }
    else
    {
      if ( elemL_[k] != interface.elemR(id) || partL_[k] != interface.partR(id) )
      {
        interface.print();
        print_inline( facet.indices );
        printf("id = %d, elemL = %lu, partL = %lu\n",id,elemL_[k],partL_[k]);
      }
      avro_assert( elemL_[k] == interface.elemR(id) && partL_[k] == interface.partR(id) );
      elemR_[k] = elemL;
      partR_[k] = partL;
    }
  }
}


template<typename type>
int
PartitionBoundary<type>::find( const ElementIndices& f ) const
{
  std::map<ElementIndices,index_t>::const_iterator it = facet_.find(f);
  if (it == facet_.end()) return -1;
  return int(it->second);
}

#define TAG_PARTITION_BND_DATA 501
#define TAG_PARTITION_BND_ELEML 502
#define TAG_PARTITION_BND_PARTL 503
#define TAG_PARTITION_BND_ELEMR 504
#define TAG_PARTITION_BND_PARTR 505

template<typename type>
void
PartitionBoundary<type>::send( index_t receiver ) const
{
  std::vector<index_t> data;
  for (std::map<ElementIndices,index_t>::const_iterator it=facet_.begin();it!=facet_.end();++it)
  {
    const ElementIndices& f = it->first;
    for (index_t j=0;j<f.indices.size();j++)
      data.push_back( f.indices[j] );
  }

  mpi::send( mpi::blocking{} , data , receiver , TAG_PARTITION_BND_DATA );
  mpi::send( mpi::blocking{} , elemL_ , receiver , TAG_PARTITION_BND_ELEML );
  mpi::send( mpi::blocking{} , partL_ , receiver , TAG_PARTITION_BND_PARTL );
  mpi::send( mpi::blocking{} , elemR_ , receiver , TAG_PARTITION_BND_ELEMR );
  mpi::send( mpi::blocking{} , partR_ , receiver , TAG_PARTITION_BND_PARTR );
}

template<typename type>
void
PartitionBoundary<type>::receive( index_t sender )
{

  std::vector<index_t> F = mpi::receive<std::vector<index_t>>( sender , TAG_PARTITION_BND_DATA );
  elemL_ = mpi::receive<std::vector<index_t>>( sender , TAG_PARTITION_BND_ELEML );
  partL_ = mpi::receive<std::vector<index_t>>( sender , TAG_PARTITION_BND_PARTL );
  elemR_ = mpi::receive<std::vector<index_t>>( sender , TAG_PARTITION_BND_ELEMR );
  partR_ = mpi::receive<std::vector<index_t>>( sender , TAG_PARTITION_BND_PARTR );

  index_t nb_facet = F.size()/(this->number()+1);
  avro_assert( nb_facet == elemL_.size() );
  for (index_t k=0;k<nb_facet;k++)
  {
    ElementIndices facet;
    facet.dim = this->number();
    facet.indices.assign( &F[(this->number()+1)*k] , &F[(this->number()+1)*(k+1)] );
    facet_.insert( {facet,facet_.size()} );
  }
}

template<typename type>
void
PartitionBoundary<type>::print() const
{
  for (std::map<ElementIndices,index_t>::const_iterator it=facet_.begin();it!=facet_.end();++it)
  {
    const ElementIndices& f = it->first;
    index_t k = it->second;
    avro_assert( k < elemL_.size() && k < partL_.size() && k < elemR_.size() & k < partR_.size() );
    print_inline( f.indices , "facet with elemL = " + std::to_string(elemL_[k]) + ", partL = " +
                                                      std::to_string(partL_[k]) + ", elemR = " +
                                                      std::to_string(elemR_[k]) + ", partR = " +
                                                      std::to_string(partR_[k]) + ": " );
  }
}

template class Partition<Simplex>;
template class Topology_Partition<Simplex>;
template class PartitionBoundary<Simplex>;

} // avro

#endif
