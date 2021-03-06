//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/mpi.hpp"
#include "common/process.h"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parallel.h"

#include "element/simplex.h"

#include "geometry/entity.h"

#include "graphics/application.h"

#include "library/meshb.h"
#include "library/metric.h"
#include "library/spacetime.h"

#include "mesh/boundary.h"
#include "mesh/facets.h"
#include "mesh/mpi_tags.h"
#include "mesh/partition.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include "numerics/mat.h"

#include <fstream>
#include <list>
#include <unordered_map>
#include <time.h>
#include <unistd.h>

#if AVRO_MPI
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
#endif

namespace avro
{

#if AVRO_MPI

#define TAG_METRIC_INDICES 301
#define TAG_METRIC_ENTRIES 302

#define TAG_FIXED_GLOBALS 401
#define TAG_GLOBAL_KEYS 402
#define TAG_GLOBAL_VALUES 403
#define TAG_POINT_OFFSET 404
#define TAG_BUDDY 601

void
pprintf(const char* format , ... ) {
  va_list args;
  va_start(args,format);
  index_t rank = mpi::rank();
  printf("[processor %2lu]: ",rank);
  printf(format,args);
  va_end(args);
}

void
avro_msgp( const std::string& msg ) {
  index_t rank = mpi::rank();
  if (rank == 0) {
    printf("[processor %3lu]: %s\n",index_t(0),msg.c_str());
    for (index_t k = 1; k < index_t(mpi::size()); k++) {
      std::vector<char> msgc = mpi::receive<std::vector<char>>( k , TAG_MISC );
      printf("[processor %3lu]: %s\n",k,msgc.data());
    }
    fflush(stdout);
  }
  else {
    std::vector<char> msgc(msg.begin(),msg.end());
    msgc.push_back('\0');
    mpi::send( mpi::blocking{} , msgc , 0 , TAG_MISC );
  }
  mpi::barrier();
}

// useful for debugging for problems in [0,1]^d
inline bool
check_geometry( const Points& points ) {
  real_t tol = 1e-12;
  index_t nb_error = 0;
  for (index_t k = 0; k < points.nb(); k++) {
    bool geometry = false;
    for (coord_t d = 0; d < points.dim(); d++) {
      if (fabs(points[k][d]-1.0) < tol) geometry = true;
      if (fabs(points[k][d]    ) < tol) geometry = true;
    }

    if (geometry && points.entity(k) == nullptr ) {
      printf("point %lu is on geometry but has no entity!\n",k);
      points.print(k,true);
      nb_error++;
    }
    else if (!geometry && points.entity(k) != nullptr) {
      printf("point %lu is interior but has an entity!\n",k);
      points.print(k,true);
      nb_error++;
    }
  }
  return nb_error == 0;
}

template<typename type>
AdaptationManager<type>::AdaptationManager( const Topology<type>& topology , const std::vector<VertexMetric>& metrics , AdaptationParameters& params ) :
  params_(params),
  topology_( topology.points().dim() , topology.points().udim() , topology.number() ),
  analytic_(nullptr) {

  // save some mpi stuff
  rank_ = mpi::rank();

  // initialize our working topology (may require partitioning)
  initialize(topology,metrics);
}

template<typename type>
void
AdaptationManager<type>::set_analytic( library::MetricField_Analytic* analytic ) {
  analytic_ = analytic;
}

template<typename type>
void
AdaptationManager<type>::send_metrics( index_t receiver , const std::vector<index_t>& indices , const std::vector<VertexMetric>& metrics , bool global_flag ) {

  // number of unique entries in a metric
  index_t nb_entry = topology_.number()*(topology_.number()+1)/2;
  std::vector<index_t> global_indices( indices.size() , 0 );
  std::vector<real_t> metric_data( indices.size()*nb_entry , 0.0 );

  index_t entry = 0;
  if (global_flag) {
    // indices are global, we we need to retrieve the local metrics
    for (index_t k = 0; k < indices.size(); k++) {
      index_t g = indices[k];
      const VertexMetric& m = metrics[g];
      for (int i = 0; i < m.n(); i++)
      for (int j = i; j < m.n(); j++)
        metric_data[entry++] = m(i,j);
    }
    global_indices = indices;
  }
  else {
    // indices are local, so we need to convert to global
    for (index_t k=0;k<indices.size();k++) {
      const VertexMetric& m = metrics[indices[k]];
      for (int i = 0; i < m.n(); i++)
      for (int j = i; j < m.n(); j++)
        metric_data[entry++] = m(i,j);
      global_indices[k] = topology_.points().global( indices[k] );
    }
  }

  // send the data
  mpi::send( mpi::blocking{} , global_indices , receiver , TAG_METRIC_INDICES );
  mpi::send( mpi::blocking{} , metric_data , receiver , TAG_METRIC_ENTRIES );
}

template<typename type>
void
AdaptationManager<type>::receive_metrics( index_t sender , bool overwrite ) {

  // number of unique entries in a metric
  index_t n = topology_.number()*(topology_.number()+1)/2;

  // option to overwrite means all metrics will be cleared
  if (overwrite)
    metrics_.resize( topology_.points().nb() , VertexMetric(topology_.number()) );

  // receive the data
  std::vector<index_t> global_indices = mpi::receive< std::vector<index_t> >( sender , TAG_METRIC_INDICES );
  std::vector<real_t> metric_data = mpi::receive< std::vector<real_t> >( sender , TAG_METRIC_ENTRIES );

  // convert the global_indices to a map so we can look things up faster
  std::map<index_t,index_t> global;
  for (index_t k=0;k<global_indices.size();k++)
    global.insert( {global_indices[k] , k} );

  for (index_t k=0;k<topology_.points().nb();k++)
  {
    std::map<index_t,index_t>::iterator it = global.find( topology_.points().global(k) );
    avro_assert( it != global.end() );

    // determine which index this is in the list
    index_t idx = it->second;

    // unpack the metric
    VertexMetric m( topology_.number() );
    index_t entry = 0;
    for (int i=0;i<m.n();i++)
    for (int j=i;j<m.n();j++)
      m(i,j) = metric_data[ idx*n + entry++ ];
    metrics_[k] = m;
  }
}

template<typename type>
void
AdaptationManager<type>::initialize(const Topology<type>& topology , const std::vector<VertexMetric>& metrics ) {

  index_t nb_rank = mpi::size();
  index_t nb_partition = nb_rank;

  if (!params_["partitioned"]) {
    // extract the entities so we can assign them to the partitions
    std::vector<Entity*> entities;
    for (index_t k = 0; k < topology.points().nb(); k++) {
      if (topology.points().entity(k)==nullptr) continue;
      entities.push_back( topology.points().entity(k) );
    }
    uniquify(entities);

    if (rank_ == 0) {

      // partition the input topology
      Partition<type> partition(topology);
      partition.compute( nb_partition );

      // extract the partitions
      std::vector<std::shared_ptr<Topology_Partition<type>>> pieces(nb_partition);
      partition.get(pieces);

      // steal partition 0 with the metrics into our working partition
      topology_.TopologyBase::copy( *pieces[0] );
      pieces[0]->points().copy( topology_.points() );
      pieces[0]->set_entities( entities );
      metrics_.resize( pieces[0]->points().nb() );
      for (index_t j = 0; j < pieces[0]->points().nb(); j++)
        metrics_[j] = metrics[ pieces[0]->points().global(j) ];

      // send off the remaining partitions to the processors
      for (index_t k = 1; k < nb_rank; k++) {
        pieces[k]->set_entities(entities);
        pieces[k]->send(k);

        std::vector<index_t> indices(pieces[k]->points().nb());
        for (index_t j = 0; j < pieces[k]->points().nb(); j++)
          indices[j] = pieces[k]->points().global(j);
        send_metrics( k , indices , metrics , true );
      }
    }
    else {
      // receive our partition
      topology_.set_entities(entities);
      topology_.receive(0);

      // receive our set of metrics
      receive_metrics(0,true);
    }
    topology_.set_entities(entities);
  }
  else {

    // copy the input topology into the working topology
    topology_.TopologyBase::copy( topology );
    topology.points().copy( topology_.points() );

    // copy the input metrics
    metrics_ = metrics;
  }
  mpi::barrier();
}

template<typename type>
class PartitionGraph {
public:
  PartitionGraph( Topology<type>& partition , const PartitionBoundary<type>& boundary , const std::vector<index_t>& vtxdist , const std::vector<VertexMetric>& metrics ) :
    partition_(partition),
    boundary_(boundary),
    vtxdist_(vtxdist) {
    avro_assert( vtxdist_.size() == index_t(mpi::size()+1) );
    build(metrics);
  }

  void build( const std::vector<VertexMetric>& metrics ) {

    index_t rank = mpi::rank();

    edges_.clear();
    adjncy_.clear();
    xadj_.clear();
    egwt_.clear();
    quality_.clear();

    mpi::barrier();

    // determine which elements may have been fixed
    std::set<index_t> fixed;
    for (index_t k = 0; k < partition_.points().nb(); k++) {

      if (!partition_.points().fixed(k)) continue;

      std::vector<index_t> ball;
      partition_.inverse().ball( k , ball );

      for (index_t j = 0; j < ball.size(); j++)
        fixed.insert( global_elem( rank , ball[j] ) );
    }
    //avro_msgp( "number of fixed elems = " + std::to_string(fixed.size()) );

    std::vector<index_t> elem_age( partition_.nb() , 1 );
    for (index_t k = 0; k < partition_.nb(); k++) {
      for (index_t j = 0; j < partition_.nv(k); j++)
        elem_age[k] += partition_.points().age( partition_(k,j) );
    }
    index_t max_age0 = *std::max_element( elem_age.begin() , elem_age.end() );

    index_t max_age = (index_t) mpi::all_reduce( int(max_age0) , mpi::max{} );

    index_t penalty = 10; // heuristic
    for (index_t k = 0; k < elem_age.size(); k++)
	    elem_age[k] = index_t( 10.0*elem_age[k]/max_age ) + 1;

    // extract the adjacencies from the partition
    std::map<index_t,index_t> local;
    for (index_t k = 0; k < partition_.nb(); k++) {
      for (index_t j = 0; j < partition_.neighbours().nfacets(); j++) {

        int n = partition_.neighbours()(k,j);
        if (n < 0) continue;      // do not create adjacency information for boundary elements
        if (n < int(k)) continue; // only create an edge once

        index_t g0 = global_elem( rank , k );
        index_t g1 = global_elem( rank , index_t(n) );

        edges_.push_back( g0 );
        edges_.push_back( g1 );
      }
    }

    // create the global to local information
    for (index_t k = 0; k < partition_.nb(); k++)
      local.insert( {global_elem(rank,k),k} );

    // add the edges corresponding to the boundary facets
    std::set<index_t> ghost;
    for (index_t k = 0; k < boundary_.nb(); k++) {

      index_t g0 = global_elem( boundary_.partL(k) , boundary_.elemL(k) );
      index_t g1 = global_elem( boundary_.partR(k) , boundary_.elemR(k) );

      if (boundary_.partL(k) == rank ) {

        avro_assert( boundary_.partR(k) != rank );
        avro_assert( local.find(g0) != local.end() );
        edges_.push_back( g0 );
        edges_.push_back( g1 );

        // keep track of ghost elements
        ghost.insert( g1 );
      }
      else {

        avro_assert( boundary_.partR(k) == rank );
        avro_assert( local.find(g1) != local.end() );
        edges_.push_back( g1 );
        edges_.push_back( g0 );

        // keep track of ghost elements
        ghost.insert(g0);
      }
    }
    nb_vert_ = partition_.nb(); // number of vertices local to this processor

    // assign weights on the node of the graph (elements in the mesh)
    // (these may not actually get used by the graph partitioner)
    assign_weights(metrics);

    // convert to csr
    std::vector< std::vector<index_t> > node2node(nb_vert_);
    std::vector< std::vector<index_t> > node2node_weight( nb_vert_ );
    for (index_t k = 0; k < edges_.size()/2; k++) {

      index_t elem0 = edges_[2*k  ];
      index_t elem1 = edges_[2*k+1];
      avro_assert( elem0 != elem1 );

      // determine the weight on this edge from the age of the elements, or if they are fixed
      #if 1
      index_t q0 = (ghost.find(elem0) == ghost.end()) ? elem_age[local.at(elem0)] : penalty;
      index_t q1 = (ghost.find(elem1) == ghost.end()) ? elem_age[local.at(elem1)] : penalty;
      #else
      index_t q0 = (ghost.find(elem0) == ghost.end()) ? 1 : penalty;
      index_t q1 = (ghost.find(elem1) == ghost.end()) ? 1 : penalty;
      #endif

      bool f0 = fixed.find(elem0) != fixed.end() || ghost.find(elem0) != ghost.end();
      bool f1 = fixed.find(elem1) != fixed.end() || ghost.find(elem1) != ghost.end();
      if (f0) q0 = penalty;
      if (f1) q1 = penalty;

      // only add adjacency information for vertices (elements) on this processor
      if (ghost.find(elem0) == ghost.end()) {
        node2node[local.at(elem0)].push_back(elem1);
        node2node_weight[local.at(elem0)].push_back( q0 );
      }
      if (ghost.find(elem1) == ghost.end()) {
        node2node[local.at(elem1)].push_back(elem0);
        node2node_weight[local.at(elem1)].push_back( q1 );
      }
    }

    // convert to the parmetis arrays
    xadj_.resize( nb_vert_+1 );
    xadj_[0] = 0;
    adjncy_.clear();
    egwt_.clear();
    for (index_t k = 0; k < nb_vert_; k++) {

      // set the adjacency pointers
      xadj_[k+1] = xadj_[k] + node2node[k].size();

      std::vector<index_t> adjlist( node2node[k].begin(),node2node[k].end() );
      std::vector<index_t> wgtlist( node2node_weight[k].begin() , node2node_weight[k].end() );

      avro_assert( adjlist.size() == wgtlist.size() );
      for (index_t i = 0; i < adjlist.size(); i++) {
        adjncy_.push_back(adjlist[i]);
        egwt_.push_back(wgtlist[i]);
      }
    }
    mpi::barrier();

  }

  void check() {

    // this function should only be called for small problems
    index_t N = vtxdist_[mpi::size()];
    avro_assert_msg( N < 1000 , "this should only be used for debugging on small problems");
    if (mpi::rank() == 0) {

      matd<real_t> graph(N,N);

      index_t n = 0;
      for (index_t k = 0; k < xadj_.size()-1; k++) {
        index_t p = xadj_[k];
        index_t q = xadj_[k+1];
        for (index_t j = p; j < q; j++)
          graph(n,adjncy_[j]) = 1.0;
        n++;
      }

      for (index_t r = 1; r < mpi::size(); r++) {

        std::vector<index_t> xadj = mpi::receive<std::vector<index_t>>(r,TAG_CELL_FIRST);
        std::vector<index_t> adjncy = mpi::receive<std::vector<index_t>>(r,TAG_CELL_LAST);
        for (index_t k = 0; k < xadj.size()-1; k++) {
          index_t p = xadj[k];
          index_t q = xadj[k+1];
          for (index_t j = p; j < q; j++)
            graph(n,adjncy[j]) = 1.0;
          n++;
        }
      }

      index_t nb_error = 0;
      for (index_t k = 0; k < graph.m(); k++)
      for (index_t i = 0; i < graph.n(); i++) {
        if (graph(k,i) != graph(i,k)) {
          printf( "non-symmetric adjacency matrix in entry (%lu,%lu) = %g but (%lu,%lu) = %g\n" ,k,i , graph(k,i) , i,k,graph(i,k) );
          nb_error++;
        }
      }
      avro_assert_msg( nb_error == 0 , "non-symmetric adjacency matrix");
    }
    else {
      mpi::send( mpi::blocking{} , xadj_ , 0 , TAG_CELL_FIRST );
      mpi::send( mpi::blocking{} , adjncy_ , 0 , TAG_CELL_LAST );
    }
    mpi::barrier();
  }

  void assign_weights( const std::vector<VertexMetric>& metrics ) {

    real_t alpha = 10.0;
    real_t beta = 0.1;

    if (metrics.size() == 0) {
      // it's possible to have partitions with no points if there are too many partitions
      vgwt_.resize( nb_vert_ , 0. );
      return;
    }

    // create a metric field to make computations easier
    MetricAttachment attachment( partition_.points() , metrics );
    MetricField<type> metric( partition_ , attachment );

    // assign a weight based on digonnet 2019
    vgwt_.resize( nb_vert_ , 0.0 );
    quality_.resize( nb_vert_ );
    for (index_t k = 0; k < partition_.nb(); k++) {
      real_t q = metric.quality(partition_,k);
      quality_[k] = q;
      avro_assert( q >= 0.0 );
      vgwt_[k] = alpha*( 1. + beta*fabs( 1. - 1./q ) );
    }

  }

  index_t global_elem( index_t partition , index_t elem ) const {
    avro_assert_msg( partition < vtxdist_.size()-1 , "partition = %lu" , partition );
    return vtxdist_[ partition ] + elem;
  }

  void repartition( std::vector<index_t>& repartition , index_t nb_part = 0 ) {

    // setup the parmetis version of the adjacency graph
    std::vector<PARM_INT> vtxdist(vtxdist_.begin(),vtxdist_.end());
    std::vector<PARM_INT> xadj(xadj_.begin(),xadj_.end());
    std::vector<PARM_INT> adjncy(adjncy_.begin(),adjncy_.end());
    std::vector<PARM_INT> vwgt( vgwt_.begin() , vgwt_.end() );
    std::vector<PARM_INT> ewgt( egwt_.begin() , egwt_.end() );
    PARM_INT* pvwgt = NULL; // vwgt.data();
    PARM_INT *padjwgt = ewgt.data();// NULL;
    PARM_INT wgtflag = 1; // 0 for no weights, 1 for edge weights, 2 for vertex weights, 3 for both
    PARM_INT edgecut = 0;
    PARM_INT bias = 0;
    PARM_INT ncon = 1;
    PARM_INT nparts = nb_part;
    if (nb_part == 0) nparts = mpi::size();

    avro_assert( ewgt.size() == adjncy.size() );

    std::vector<PARM_REAL> tpwgts(ncon*nparts,1./nparts);
    std::vector<PARM_REAL> ubvec(ncon,1.001);//1.05); // 1.001 if segfault
    PARM_INT options[4];
    options[0] = 1;
    options[1] = 0;
    options[2] = 42;
    //options[3] = PARMETIS_PSR_COUPLED; // only used for AdaptiveRepart

    // partition the graph!
    std::vector<PARM_INT> part(nb_vert_,mpi::rank());
    MPI_Comm comm = MPI_COMM_WORLD;
    #if 1
    int result = ParMETIS_V3_PartKway( vtxdist.data() , xadj.data() , adjncy.data() ,
                            pvwgt , padjwgt, &wgtflag,
                            &bias , &ncon , &nparts,
                            tpwgts.data() , ubvec.data(),
                            options,
                            &edgecut,
                            part.data(),
                            &comm);
    #elif 0 // this only works for d <= 3

    PARM_INT dim = partition_.points().dim();
    std::vector<PARM_REAL> coordinates( partition_.nb() * dim );

    for (index_t k = 0; k < partition_.nb(); k++) {
      index_t nv = partition_.nv(k);
      for (index_t j = 0; j < nv; j++)
      for (coord_t d = 0; d < dim; d++)
        coordinates[dim*k+d] += partition_.points()[ partition_(k,j) ][d] / nv;
    }

    int result = ParMETIS_V3_PartGeomKway( vtxdist.data() , xadj.data() , adjncy.data() ,
                            pvwgt , padjwgt, &wgtflag,
                            &bias , &dim , coordinates.data() , &ncon , &nparts,
                            tpwgts.data() , ubvec.data(),
                            options,
                            &edgecut,
                            part.data(),
                            &comm);

    #else
    PARM_REAL itr = 1000000.;
    index_t mem_vertex = (partition_.number()+1)*sizeof(index_t);
    std::vector<PARM_INT> vsize(nb_vert_,mem_vertex);
    int result = ParMETIS_V3_AdaptiveRepart( vtxdist.data() , xadj.data() , adjncy.data() ,
                            vwgt.data(), vsize.data(), padjwgt, &wgtflag,
                            &bias , &ncon , &nparts,
                            tpwgts.data() , ubvec.data(), &itr,
                            options,
                            &edgecut,
                            part.data(),
                            &comm);
    #endif
    avro_assert( result == METIS_OK );
    if (mpi::rank() == 0)
      printf("partitioned mesh into %lu domains, final edge cut = %d\n",nb_part,edgecut);

    // store the partition result
    repartition.assign( part.begin() , part.end() );
    mpi::barrier();
  }

private:
  Topology<type>& partition_;
  const PartitionBoundary<type>& boundary_;

  std::vector<index_t> edges_;
  std::vector<index_t> adjncy_;
  std::vector<index_t> xadj_;

  std::vector<real_t> vgwt_;
  std::vector<real_t> egwt_;

  const std::vector<index_t>& vtxdist_;

  std::vector<real_t> quality_;

  index_t nb_vert_;
};

class ProcessorGraph {
public:
  ProcessorGraph( const std::vector<index_t>& element_offset , const std::set< std::pair<index_t,index_t> >& pairs ) :
    element_offset_(element_offset),
    previous_pairs_(pairs)
  {}

  template<typename type>
  void match( const PartitionBoundary<type>& interface ) {

    avro_assert( mpi::size()%2 == 0 );

    std::set< std::pair<index_t,index_t> > E;
    const std::map<ElementIndices,index_t>& facets = interface.facets();

    std::map<ElementIndices,index_t>::const_iterator it;
    for (it = facets.begin(); it != facets.end(); ++it) {

      index_t k = it->second;

      index_t partL = interface.partL(k);
      index_t partR = interface.partR(k);
      avro_assert( partL != partR );

      // always keep processors with lower id's to the left
      if (partL > partR) std::swap(partL,partR);
      E.insert( {partL,partR} );
    }

    #if 0
    PerfectMatching graph( mpi::size() , E.size() );
    graph.options.verbose = false;
    std::vector<index_t> edges;
    std::vector<index_t> edge_ids;
    for (std::set< std::pair<index_t,index_t> >::iterator it=E.begin();it!=E.end();++it)
    {
      index_t p0 = it->first;
      index_t p1 = it->second;
      edges.push_back(p0);
      edges.push_back(p1);

      // TODO determine a better cost function
      // likely the number of facets on this interface
      index_t nb0 = element_offset_[p0+1] - element_offset_[p0];
      index_t nb1 = element_offset_[p1+1] - element_offset_[p1];
      int cost = ( nb0 + nb1 );

      // heavily penalize edges in the matching that would repeat a previous matching
      if (previous_pairs_.find( {p0,p1} ) != previous_pairs_.end())
        cost = 100*cost;

      edge_ids.push_back( graph.AddEdge( p0 , p1 , cost ) );
    }

    // solve the graph matching problem
    graph.Solve();
    for (index_t k=0;k<edge_ids.size();k++)
    {
      index_t id = edge_ids[k];
      if (graph.GetSolution(id) == 0) continue; // not an edge in the matching
      pairs_.push_back( {edges[2*k],edges[2*k+1]} );
    }
    #else
    printf("blossom5 was removed since this was no longer needed");
    #endif
  }

  index_t nb_pairs() const { return pairs_.size(); }
  std::vector<std::pair<index_t,index_t>> pairs() const { return pairs_; }

  std::pair<index_t,index_t> get_pair( index_t k ) { return pairs_[k]; }

private:
  const std::vector<index_t>& element_offset_;
  std::vector< std::pair<index_t,index_t> > pairs_;
  const std::set<std::pair<index_t,index_t>>& previous_pairs_;
};

template<typename type>
class MigrationChunk : public Topology_Partition<type>
{
public:
  MigrationChunk( coord_t dim , coord_t udim , coord_t number ) :
    Topology_Partition<type>(dim,udim,number)
  {}

  void extract( const Topology_Partition<type>& topology , const std::vector<index_t>& partition , index_t rank )
  {
    // first add any element that stays in our partition
    // keep track of the global indices that have been added
    partition_.reserve( topology.nb() );
    Topology<type>::reserve( topology.nb() );
    for (index_t k=0;k<topology.nb();k++)
    {
      if (partition[k] != rank) continue;

      std::vector<index_t> s = topology.get(k); // in local topology indexing
      for (index_t j=0;j<s.size();j++)
      {
        index_t p = topology.points().global(s[j]);
        if (global_.find(p) == global_.end())
        {
          // create the point
          index_t q = points_.nb();
          index_t m = s[j];
          points_.create( topology.points()[m] );
          points_.set_entity( q , topology.points().entity(m) );
          points_.set_param( q , topology.points().u(m) );
          points_.set_global( q , topology.points().global(m) );
          points_.set_fixed( q , topology.points().fixed(m) );
          global_.insert( {p,q} );
        }
        s[j] = global_.at(p);
      }
      this->add(s.data(),s.size());
      partition_.push_back(rank);
    }
    partition_.shrink_to_fit();
    Topology<type>::shrink_to_fit();

    avro_assert_msg( global_.size() == points_.nb() , "|global| = %lu, |points| = %lu" , global_.size() , points_.nb() );
  }

  void add_to( Topology<type>& topology ) {
    std::map<index_t,index_t> global;
    for (index_t k = 0; k < topology.points().nb(); k++)
      global.insert( {topology.points().global(k) , k} );

    topology.points().reserve( topology.points().nb() + points_.nb() );
    metric_.reserve( metric_.size() + points_.nb() );
    for (index_t k = 0; k < this->nb(); k++) {

      std::vector<index_t> s = this->get(k);
      for (index_t j = 0; j < s.size(); j++) {

        // check if this index already exists in the global list
        index_t p = points_.global(s[j]);
        index_t q;
        if (global.find(p) == global.end())  {

          // this point is added to the topology, so the metric will be needed
          metric_.push_back( p );

          // create the point
          q = topology.points().nb();
          index_t m = s[j];
          topology.points().create( points_[m] );
          topology.points().set_entity( q , points_.entity(m) );
          topology.points().set_param( q , points_.u(m) );
          topology.points().set_global( q , points_.global(m) );
          topology.points().set_fixed( q , points_.fixed(m) );
          global.insert( {p,q} );
        }
        else
          q = global.at(p);
        s[j] = q;
      }
      topology.add( s.data() , s.size() );
    }
    topology.points().shrink_to_fit();
    metric_.shrink_to_fit();
  }

  const std::vector<index_t>& partition() const { return partition_; }
  const std::vector<index_t>& metric() const { return metric_; }

private:
  using Topology_Partition<type>::points_;
  std::unordered_map<index_t,index_t> global_;
  std::vector<index_t> partition_;
  std::vector<index_t> metric_;
};

template<typename type>
void
AdaptationManager<type>::exchange( std::vector<index_t>& repartition ) {

  index_t nb_rank = mpi::size();

  if (rank_ == 0) printf("--> exchanging elements\n");
  mpi::barrier();

  // determine if there will be zero elements remaining on this processor
  index_t nb_keep = 0;
  for (index_t k = 0; k < repartition.size(); k++) {
    if (repartition[k] == rank_) nb_keep++;
  }
  if (nb_keep == 0) {

    // find the single best element to keep
    MetricAttachment attachment( topology_.points() , metrics_ );
    MetricField<type> field(topology_,attachment);

    index_t kmax = 0;
    real_t qmax = 0;
    for (index_t k = 0; k < topology_.nb(); k++) {
      real_t q = field.quality( topology_, k );
      if (q > qmax) {
        qmax = q;
        kmax = k;
      }
    }
    repartition[kmax] = rank_;
  }

  // print a message with all the send information
  std::vector<index_t> nb_send( nb_rank , 0 );
  for (index_t k = 0; k < repartition.size(); k++)
    nb_send[repartition[k]]++;
  std::string msg = "exchange: ( ";
  for (index_t j = 0; j < nb_send.size(); j++)
    msg += std::to_string(nb_send[j]) + " ";
  msg += ")";
  avro_msgp( msg );

  coord_t dim = topology_.points().dim();
  coord_t udim = topology_.points().udim();
  coord_t number = topology_.number();

  clock_t TIME0, TIME1;

  // accumulate the elements to send away
  TIME0 = clock();
  MigrationChunk<type> goodbye(dim,udim,number);
  goodbye.set_entities(topology_.entities());
  for (index_t j = 0; j < nb_rank; j++) {
    if (j == rank_) continue;
    goodbye.extract( topology_ , repartition , j );
  }
  mpi::barrier();
  TIME1 = clock();
  if (rank_ == 0) printf("--> extract time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // receive chunks on the root processor so we can send the combined chunks away
  TIME0 = clock();
  std::vector< std::shared_ptr<MigrationChunk<type>> > pieces( nb_rank );
  if (rank_ == 0) {

    // initialize the pieces
    for (index_t k = 0; k < nb_rank; k++) {
      pieces[k] = std::make_shared<MigrationChunk<type>>(dim,udim,number);
      pieces[k]->set_entities(topology_.entities());
    }

    // extract the pieces from the root processor
    for (index_t j = 1; j < nb_rank; j++)
      pieces[j]->extract( goodbye , goodbye.partition() , j );

    for (index_t k = 1; k < nb_rank; k++) {

      // receive the chunk from processor k
      MigrationChunk<type> chunk(dim,udim,number);
      chunk.set_entities( topology_.entities() );
      chunk.receive( k );
      std::vector<index_t> partition_k = mpi::receive<std::vector<index_t>>(k,TAG_MISC);
      avro_assert( partition_k.size() == chunk.nb() );

      // add the elements to the appropriate pieces
      for (index_t j = 0; j < nb_rank; j++)
        pieces[j]->extract( chunk , partition_k , j );
    }
  }
  else {
    // send the elements we don't need anymore to the root processor
    goodbye.send(0);
    mpi::send( mpi::blocking{} , goodbye.partition() , 0 , TAG_MISC );
  }
  mpi::barrier();
  TIME1 = clock();
  if (rank_ == 0) printf("--> chunk sending to root time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // send the chunks to the corresponding processors
  TIME0 = clock();
  MigrationChunk<type> chunk(dim,udim,number);
  chunk.set_entities( topology_.entities() );
  if (rank_ == 0) {

    // accumulate the pieces
    for (index_t k = 1; k < nb_rank; k++)
      pieces[k]->send(k);

    // append the root piece into the chunk
    chunk.TopologyBase::copy( *pieces[0].get() );
    pieces[0]->points().copy(chunk.points());
  }
  else {
    // receive the chunks we need to add
    chunk.receive(0);
  }
  mpi::barrier();
  TIME1 = clock();
  if (rank_ == 0) printf("--> chunk send to proc time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));


  std::vector<index_t> removals;
  for (index_t k = 0; k < repartition.size(); k++) {
    if (repartition[k] != rank_)
      removals.push_back( k );
  }

  std::unordered_map<index_t,index_t> metrics_owned;
  for (index_t k = 0; k < topology_.points().nb(); k++)
    metrics_owned.insert( {topology_.points().global(k),k} );

  // add the received chunk to the topology
  TIME0 = clock();
  index_t nb_elem = topology_.nb() - removals.size();
  topology_.reserve( topology_.nb() + chunk.nb() );
  chunk.add_to( topology_ );
  topology_.shrink_to_fit();
  TIME1 = clock();
  if (rank_ == 0) printf("--> chunk addition time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  crust_.clear();
  crust_.resize( chunk.nb() );
  for (index_t k = 0; k < chunk.nb(); k++)
    crust_[k] = nb_elem + k;
    //crust_.push_back( nb_elem + k );

  TIME0 = clock();
  topology_.remove_elements( removals );
  for (index_t k = 0; k < crust_.size(); k++)
    avro_assert_msg( crust_[k] < topology_.nb(), "crust = %lu, |topology| = %lu on rank %lu",crust_[k],topology_.nb(),rank_);
  TIME1 = clock();
  if (rank_ == 0) printf("--> element removal time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // determine which metrics need to be added
  TIME0 = clock();
  std::set<index_t> idx( chunk.metric().begin() , chunk.metric().end() );
  std::vector<index_t> metric;
  metric.reserve( chunk.metric().size() );
  for (index_t k = 0; k < chunk.metric().size(); k++) {
    index_t p = chunk.metric()[k];
    if (metrics_owned.find(p) == metrics_owned.end())
      metric.push_back(p);
  }
  metric.shrink_to_fit();

  // accumulate the requested metrics on the root processor
  std::vector<std::set<index_t>> metric_send_idx(nb_rank);
  std::vector<index_t> metric_requests;
  if (rank_ == 0)
  {
    // initialize the full list of metrics on the root processor
    metric_requests = metric;

    metric_requests.reserve( metric.size() * nb_rank );

    // append the metric requests received from other processors
    for (index_t k = 1; k < nb_rank; k++) {
      std::vector<index_t> metric_k = mpi::receive<std::vector<index_t>>( k , TAG_MISC );
      std::copy( metric_k.begin() , metric_k.end() , std::inserter(metric_send_idx[k],metric_send_idx[k].end() ) );
      for (index_t j = 0; j < metric_k.size(); j++)
        metric_requests.push_back(metric_k[j]);
    }
    metric_requests.shrink_to_fit();
    uniquify(metric_requests);
  }
  else {
    mpi::send( mpi::blocking{} , metric , 0 , TAG_MISC );
  }
  mpi::barrier();
  TIME1 = clock();
  if (rank_ == 0) printf("--> metric accumulation (on root) time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  TIME0 = clock();
  // now gather the data from all the processors
  if (rank_ == 0) {
    for (index_t k = 1; k < nb_rank; k++)
      mpi::send( mpi::blocking{} , metric_requests , k , TAG_MISC );
  }
  else {
    metric_requests = mpi::receive<std::vector<index_t>>(0,TAG_MISC);
  }
  mpi::barrier();
  TIME1 = clock();
  if (rank_ == 0) printf("--> accumulated metric receive time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // create the metric data to send for the points we own
  TIME0 = clock();
  std::vector<real_t> metric_data;
  std::vector<index_t> metric_idx;
  metric_data.reserve( metric_requests.size() * number * (number+1 ) );
  metric_idx.reserve( metric_requests.size() );
  for (index_t k = 0; k < metric_requests.size(); k++) {

    // determine if we own this metric
    index_t p = metric_requests[k];
    if (metrics_owned.find(p) == metrics_owned.end()) continue;

    index_t q = metrics_owned.at(p);
    metric_idx.push_back( metric_requests[k] );
    for (index_t j = 0; j < number; j++)
    for (index_t i = j; i < number; i++)
      metric_data.push_back( metrics_[q](i,j) );
  }
  metric_data.shrink_to_fit();
  metric_idx.shrink_to_fit();

  if (rank_ == 0) {

    // receive all the metric data and append it to the list (initialized by the root)
    for (index_t k = 1; k < nb_rank; k++) {

      std::vector<index_t> metric_idx_k = mpi::receive<std::vector<index_t>>( k , TAG_METRIC_INDICES );
      std::vector<real_t> metric_data_k = mpi::receive<std::vector<real_t>>( k , TAG_METRIC_ENTRIES );

      // append the data
      index_t count = 0;
      for (index_t m = 0; m < metric_idx_k.size(); m++) {
        metric_idx.push_back(metric_idx_k[m]);
        for (index_t j = 0; j < number; j++)
        for (index_t i = j; i < number; i++)
          metric_data.push_back( metric_data_k[count++] );
      }
    }
  }
  else {
    // send our metric data to the root
    mpi::send( mpi::blocking{} , metric_idx , 0 , TAG_METRIC_INDICES );
    mpi::send( mpi::blocking{} , metric_data , 0 , TAG_METRIC_ENTRIES );
  }
  mpi::barrier();

  // now receive all the data from the host
  if (rank_ == 0) {
    for (index_t k = 1; k < nb_rank; k++) {
      mpi::send( mpi::blocking{} , metric_idx , k , TAG_METRIC_INDICES );
      mpi::send( mpi::blocking{} , metric_data , k , TAG_METRIC_ENTRIES );
    }
  }
  else {
    metric_idx  = mpi::receive<std::vector<index_t>>(0,TAG_METRIC_INDICES);
    metric_data = mpi::receive<std::vector<real_t>>(0,TAG_METRIC_ENTRIES);
  }
  mpi::barrier();

  // find the data we need and append it into our array
  metrics_.reserve( metric_idx.size() + metrics_.size() );
  index_t N = number*(number+1)/2;
  for (index_t k = 0; k < metric.size(); k++) {
    bool found = false;
    for (index_t m = 0; m < metric_idx.size(); m++) {
      if (metric_idx[m] == metric[k]) {
        VertexMetric mk(number);
        index_t count = 0;
        for (index_t j = 0; j < number; j++)
        for (index_t i = j; i < number; i++)
          mk(i,j) = metric_data[ m*N + count++ ];
        metrics_.push_back(mk);
        found = true;
        break;
      }
    }
    avro_assert( found );
  }
  avro_assert( metrics_.size() == topology_.points().nb() );
  metrics_.shrink_to_fit();
  TIME1 = clock();
  if (rank_ == 0) printf("--> other metric time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // remove points that are not referenced by the topology
  TIME0 = clock();
  std::vector<index_t> point_map;
  topology_.remove_unused(&point_map); // i think this is a bottleneck but it's really hard to improve
  TIME1 = clock();
  if (rank_ == 0) printf("--> topology unused removal time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  // map the metrics too
  TIME0 = clock();
  #if 0
  for (index_t k = 0; k < point_map.size(); k++)
    metrics_.erase( metrics_.begin() + point_map[k]-k );
  #elif 0
  std::list<VertexMetric> lmetrics(metrics_.begin(),metrics_.end());
  std::list<VertexMetric>::iterator it = lmetrics.begin();
  for (index_t k = 0; k < point_map.size(); k++) {
    advance(it,point_map[k]-k);
    lmetrics.erase( it );
    it = lmetrics.begin();
  }
  metrics_.assign( lmetrics.begin() , lmetrics.end() );
  #else
  std::set<index_t> spoint_map(point_map.begin(),point_map.end());
  std::vector<VertexMetric> metrics0;
  metrics0.reserve( metrics_.size() );
  for (index_t k = 0; k < metrics_.size(); k++) {
    if (spoint_map.find(k) != spoint_map.end()) continue;
    metrics0.push_back(metrics_[k]);
  }
  metrics0.shrink_to_fit();
  metrics_.assign(metrics0.begin(),metrics0.end());
  #endif
  TIME1 = clock();
  if (rank_ == 0) printf("--> metric removal time: %4.3g sec.\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));


  avro_assert( metrics_.size() == topology_.points().nb() );
  for (index_t k = 0; k < topology_.nb(); k++) {
    for (index_t j = 0; j < topology_.nv(k); j++)
      avro_assert( topology_(k,j) < topology_.points().nb() );
  }

  mpi::barrier();
}

template<typename type>
void
AdaptationManager<type>::migrate_balance( index_t nb_part ) {

  index_t nb_rank = mpi::size();
  //if (nb_rank == 1) return;

  // migrate the interface into the interior of the partitions
  // simultaneously performing a load balance
  topology_.build_structures();

  // exchange partition boundaries to all processors
  PartitionBoundary<type> boundary(topology_.number()-1);
  boundary.compute(topology_,rank_);

  PartitionBoundary<type> interface( topology_.number()-1 );
  if (rank_ == 0) {

    // receive and accumulate all boundary information
    interface.append( boundary );
    for (index_t k = 1; k < nb_rank; k++) {
      PartitionBoundary<type> boundary_k( topology_.number()-1 );
      boundary_k.receive(k);
      interface.append(boundary_k);
    }
    avro_assert( interface.complete(element_offset_) );
    printf("--> number of interface facets = %lu\n",interface.nb());
  }
  else {
    // send our partition boundary information
    boundary.send(0);
  }
  mpi::barrier();

  if (rank_ == 0 ) {
    // send the accumulated boundary information
    for (index_t k = 1; k < nb_rank; k++)
      interface.send(k);
  }
  else
  {
    // receive the accumulated boundary information
    interface.receive(0);
  }
  mpi::barrier();

  // fill in the missing elemR/partR information
  boundary.fill( interface );

  // build up the local adjacency graph, accounting for elements (graph vertices)
  // on other processors stored in the elemR/partR information
  PartitionGraph<type> graph( topology_ , boundary , element_offset_ , metrics_ );

  // call parmetis to perform the load balance
  clock_t TIMER = clock();
  std::vector<index_t> repartition;
  graph.repartition(repartition,nb_part);
  if (rank_ == 0) time_partition_ += real_t(clock() - TIMER)/real_t(CLOCKS_PER_SEC);

  // exchange the elements between partitions
  TIMER = clock();
  exchange(repartition);
  if (rank_ == 0) time_exchange_ += real_t(clock() - TIMER)/real_t(CLOCKS_PER_SEC);
}

template<typename type>
void
AdaptationManager<type>::migrate_interface() {

  index_t nb_rank = mpi::size();

  // migrate the interface into the interior of the partitions
  // simultaneously performing a load balance
  topology_.build_structures();

  // exchange partition boundaries to all processors
  PartitionBoundary<type> boundary(topology_.number()-1);
  boundary.compute(topology_,rank_);

  PartitionBoundary<type> interface( topology_.number()-1 );
  if (rank_ == 0) {

    // receive and accumulate all boundary information
    interface.append( boundary );
    for (index_t k = 1; k < nb_rank; k++) {
      PartitionBoundary<type> boundary_k( topology_.number()-1 );
      boundary_k.receive(k);
      interface.append(boundary_k);
    }
    avro_assert( interface.complete(element_offset_) );
  }
  else {
    // send our partition boundary information
    boundary.send(0);
  }
  mpi::barrier();

  if (rank_ == 0 ) {
    // send the accumulated boundary information
    for (index_t k=1;k<nb_rank;k++)
      interface.send(k);
  }
  else {
    // receive the accumulated boundary information
    interface.receive(0);
  }
  mpi::barrier();

  // fill in the missing elemR/partR information
  boundary.fill( interface );

  // setup a buddy system
  index_t buddy = mpi::size();
  if (rank_ == 0) {

    // build up the interprocessor graph
    // with vertex weights defined be current processor work
    // compute a matching of this graph
    print_inline(element_offset_);
    ProcessorGraph graph(element_offset_,pairs_);
    graph.match(interface);

    // send off the processor pairs
    for (index_t k=0;k<graph.nb_pairs();k++) {

      std::pair<index_t,index_t> p = graph.get_pair(k);

      index_t p0 = p.first;
      index_t p1 = p.second;

      mpi::send( mpi::blocking{} , p1 , p0 , TAG_BUDDY );
      if (p1!=rank_)
        mpi::send( mpi::blocking{} , p0 , p1 , TAG_BUDDY );

      if (p0 == rank_) buddy = p1;
    }

    // save the matched pairs so we don't match them again
    const std::vector< std::pair<index_t,index_t>>& pairs = graph.pairs();
    for (auto it=pairs.begin();it!=pairs.end();++it)
      pairs_.insert(*it);
  }
  else {
    // wait to receive our processor pair
    buddy = mpi::receive<index_t>(0,TAG_BUDDY);
  }
  avro_assert( buddy < index_t(mpi::size()) );
  mpi::barrier();

  std::vector<index_t> discard;
  std::set<index_t> interface_point;
  if (rank_ < buddy) {

    // loop through all interface facets, finding the ones this processor owns
    const std::map<ElementIndices,index_t>& facets = interface.facets();
    std::map<ElementIndices,index_t>::const_iterator it;
    std::vector<index_t> elems;

    for (it = facets.begin(); it != facets.end(); ++it) {

      const ElementIndices& f = it->first;
      index_t k = it->second;

      index_t partL = interface.partL(k);
      index_t partR = interface.partR(k);

      index_t elemL = interface.elemL(k);
      index_t elemR = interface.elemR(k);

      // always keep rank to the left and buddy to the right
      if (partL > partR) {
        std::swap(partL,partR);
        std::swap(elemL,elemR);
      }

      // make sure this facet is on the current processor pair
      if (partL != rank_ || partR != buddy) continue;

      // identify points which are on the interface of the processor pair
      for (index_t j = 0; j < f.indices.size(); j++)
        interface_point.insert( f.indices[j] );

      // we will discard the element on this processor
      discard.push_back(elemL);
    }
  }
  mpi::barrier();

  // list out the interface points
  std::set<index_t> pts;
  for (index_t k = 0; k < topology_.points().nb(); k++) {
    index_t p = topology_.points().global(k);
    if (interface_point.find(p) == interface_point.end()) continue;
    pts.insert( k );
  }

  // we will discard any element in the ball of an interface point
  std::vector<index_t> ball;
  for (std::set<index_t>::iterator it = pts.begin(); it != pts.end(); ++it) {
    ball.clear();
    topology_.inverse().ball( *it , ball );
    for (index_t j = 0; j < ball.size(); j++)
      discard.push_back(ball[j]);
  }
  uniquify(discard);

  // initialize all of the elements to this partition
  std::vector<index_t> repartition(topology_.nb(),rank_);

  // re-assign the discarded elements to the buddy
  for (index_t k = 0; k < discard.size(); k++)
    repartition[discard[k]] = buddy;
  mpi::barrier();

  // exchange the elements
  exchange(repartition);
}

void
get_fixed_globals( const Points& points , std::vector<index_t>& global ) {
  for (index_t k = 0; k < points.nb(); k++) {
    if (!points.fixed(k)) continue;
    global.push_back(points.global(k));
  }
}

template<typename key,typename value>
struct map_builder : public std::map<key,value> {
  map_builder( const std::vector<key>& keys , const std::vector<value>& values ) {
    avro_assert( keys.size() == values.size() );
    for (index_t k = 0; k < keys.size(); k++)
      this->insert( {keys[k] , values[k] } );
  }
};

template<typename type>
void
AdaptationManager<type>::synchronize() {

  index_t nb_rank = mpi::size();

  // synchronize the element offsets
  element_offset_.resize( nb_rank+1 , 0 );
  if (rank_ == 0) {

    std::vector<index_t> nb_elem( nb_rank , 0 );
    nb_elem[0] = topology_.nb();
    for (index_t k = 1; k < nb_rank; k++)
      nb_elem[k] = mpi::receive<index_t>( k , TAG_MISC );
    for (index_t k = 0; k < nb_rank;k++)
      element_offset_[k+1] = element_offset_[k] + nb_elem[k];
  }
  else
    mpi::send( mpi::blocking{} , topology_.nb() , 0 , TAG_MISC );
  mpi::barrier();

  if (rank_ == 0) {
    for (index_t k=1;k<nb_rank;k++)
      mpi::send( mpi::blocking{} , element_offset_ , k , TAG_MISC );
  }
  else
    element_offset_ = mpi::receive<std::vector<index_t>>( 0 , TAG_MISC );
  mpi::barrier();

  // synchronize all the global point indices between processors
  // the root processor needs to know how many vertices there are in total
  std::vector<index_t> global;
  get_fixed_globals( topology_.points() , global );
  std::map<index_t,index_t> fixed_map;
  index_t nb_points_total = 0;
  std::vector<index_t> nb_points( nb_rank , 0 );
  std::vector<index_t> pt_offset( nb_rank , 0 );
  if (rank_ == 0) {

    for (index_t j = 0; j < global.size(); j++)
      fixed_map.insert( {global[j] , j } );

    // receive the total number of vertices from all processors
    nb_points[0] = topology_.points().nb();

    nb_points_total = nb_points[0];
    std::vector<index_t> nb_interior( nb_rank , 0 );
    nb_interior[0] = nb_points[0] - global.size();
    for (index_t k = 1; k < nb_rank; k++) {

      // receive the total number of points
      nb_points[k] += index_t(mpi::receive<int>( k , TAG_MISC ));
      nb_points_total += nb_points[k];

      // receive the list of global indices
      std::vector<index_t> globals_k = mpi::receive<std::vector<index_t>>( k , TAG_FIXED_GLOBALS );
      nb_interior[k] = nb_points[k] - globals_k.size();

      for (index_t j = 0; j < globals_k.size(); j++) {
        index_t g = globals_k[j];
        if (fixed_map.find(g) != fixed_map.end()) continue;
        fixed_map.insert( {g,fixed_map.size()} );
      }
    }

    // compute the point offset for every processor
    pt_offset[0] = fixed_map.size();
    for (index_t k = 1; k < nb_rank; k++) {
      pt_offset[k] = pt_offset[k-1] + nb_interior[k-1];
    }
  }
  else {
    // send the number of points on this processor
    mpi::send( mpi::blocking{} , int(topology_.points().nb()) , 0 , TAG_MISC );

    // send list of global indices
    mpi::send( mpi::blocking{} , global , 0 , TAG_FIXED_GLOBALS );
  }
  mpi::barrier();

  if (rank_ == 0) {

    // assign the new global indices for this processor
    index_t count = fixed_map.size();
    for (index_t j = 0; j < topology_.points().nb(); j++) {
      // fixed points should still be at the beginning!
      if (!topology_.points().fixed(j)) {
        topology_.points().set_global( j , count++ );
        continue;
      }

      index_t g = topology_.points().global(j);
      topology_.points().set_global( j , fixed_map[g] );
    }

    // send global indices to all processors
    std::vector<index_t> keys,values;
    for (std::map<index_t,index_t>::const_iterator it = fixed_map.begin(); it != fixed_map.end(); ++it) {
      keys.push_back( it->first );
      values.push_back( it->second );
    }
    for (index_t k = 1; k < nb_rank; k++) {
      mpi::send( mpi::blocking{} , keys , k , TAG_GLOBAL_KEYS );
      mpi::send( mpi::blocking{} , values , k , TAG_GLOBAL_VALUES );
      mpi::send( mpi::blocking{} , pt_offset[k] , k , TAG_POINT_OFFSET );
    }
  }
  else {

    // receive the old global indices
    std::vector<index_t> keys = mpi::receive<std::vector<index_t>>(0,TAG_GLOBAL_KEYS);

    // receive new global indices
    std::vector<index_t> values = mpi::receive<std::vector<index_t>>(0,TAG_GLOBAL_VALUES);

    // receive the point offset
    index_t offset = mpi::receive<index_t>(0,TAG_POINT_OFFSET);

    // make a map!
    map_builder<index_t,index_t> global_map(keys,values);

    index_t count = 0;
    for (index_t k = 0; k < topology_.points().nb(); k++) {

      // fixed points should still be at the beginning!
      if (!topology_.points().fixed(k)) {
        topology_.points().set_global( k , offset + count++ );
        continue;
      }

      // current global index
      index_t g = topology_.points().global(k);
      avro_assert( global_map.find(g)!=global_map.end() );

      // assign the new global index
      topology_.points().set_global( k , global_map.at(g) );
    }
  }
  mpi::barrier();
}

template<typename type>
void
fix_non_manifold( Topology<type>& topology )
{
  // build up the list of elements attached to each vertex
  std::vector<std::vector<index_t>> node2elem( topology.points().nb() );
  for (index_t k=0;k<topology.nb();k++) {
    for (index_t j=0;j<topology.nv(k);j++)
      node2elem[ topology(k,j) ].push_back(k);
  }

  topology.build_structures();
  std::vector<index_t> ball;
  for (index_t k = 0; k < topology.points().nb(); k++) {
    ball.clear();
    topology.inverse().ball( k , ball );

    // if the ball doesn't match our precomputed list, then this is probably a non-manifold vertex
    if (ball.size() != node2elem[k].size() ) {
      // so we fix any vertex in the ball (TODO: can this be only the vertex k?)
      for (index_t j = 0; j < node2elem[k].size(); j++)
      for (index_t i = 0; i < topology.nv(node2elem[k].size()); i++)
        topology.points().set_fixed( topology(node2elem[k][j],i) , true );
    }
  }
}

template<typename type>
void
AdaptationManager<type>::fix_boundary() {

  // first unfix all the vertices in this partition's topology
  for (index_t k = 0; k < topology_.points().nb(); k++)
    topology_.points().set_fixed(k,false);

  // the set of points that will be fixed
  std::vector<index_t> pts;

  // fix all non-manifold vertices (i.e. vertices in which the ball does not match the actual inverse)
  fix_non_manifold(topology_);
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    if (topology_.points().fixed(k))
      pts.push_back(k);
  }
  //avro_msgp("detected "+ std::to_string(pts.size()) + " non-manifold vertices out of " + std::to_string(topology_.points().nb()) + "!");

  // compute the boundary
  Facets facets(topology_);
  facets.compute();
  std::vector<index_t> facet(topology_.number());
  index_t nb_fixed_facets = 0;
  index_t nb_bnd = 0;
  for (index_t k = 0; k < facets.nb(); k++) {

    // skip interior facets
    if (!facets.boundary(k)) continue;
    nb_bnd++;

    // skip geometry facets
    facets.retrieve(k,facet);
    Entity* entity = BoundaryUtils::geometryFacet( topology_.points() , facet.data() , facet.size());
    if (entity != nullptr) continue;
    nb_fixed_facets++;

    for (index_t j = 0; j < facet.size(); j++)
      pts.push_back(facet[j]);
  }
  uniquify(pts);
  std::sort(pts.begin(),pts.end());

  if (mpi::size() == 1 )
    avro_assert( nb_fixed_facets == 0 );

  // fix all the partition boundary points
  for (index_t j = 0; j < pts.size(); j++)
    topology_.points().set_fixed(pts[j],true);

  // move all the fixed points to the beginning of the list
  // this should also adjust the global indices
  // but we'll need to adjust the metrics
  std::map<index_t,index_t> point_map;
  topology_.move_to_front( pts , &point_map );

  // map the metrics
  coord_t number = topology_.number();
  std::vector<VertexMetric> mapped_metrics(metrics_.size(),VertexMetric(number));
  for (index_t k = 0; k < metrics_.size(); k++)
    mapped_metrics[ point_map[k] ] = metrics_[k];
  metrics_ = mapped_metrics;
  mpi::barrier();
}

template<typename type>
void
fix_crust( Topology<type>& topology , const std::vector<index_t>& crust0 ) {

  for (index_t k = 0; k < crust0.size(); k++) {
    if (crust0[k] >= topology.nb()) print_inline(crust0);
    avro_assert_msg( crust0[k] < topology.nb() , "crust elem %lu but topology has %lu elements" , crust0[k] , topology.nb() );
  }

  std::set<index_t> pts;
  for (index_t k = 0; k < topology.points().nb(); k++) {
    if (topology.points().fixed(k))
      pts.insert(k);
  }

  std::set<index_t> crust(crust0.begin(),crust0.end());
  index_t nb_fixed = 0;
  for (index_t k = 0; k < topology.nb(); k++) {
    // don't fix crust elements
    if (crust.find(k) != crust.end()) continue;

    bool fixed = true;
    for (index_t j = 0; j < topology.nv(k); j++) {
      if (pts.find(topology(k,j)) != pts.end()) {
        fixed = false;
        break;
      }
    }
    if (!fixed) continue;

    nb_fixed++;
    for (index_t j = 0; j < topology.nv(k); j++)
      topology.points().set_fixed( topology(k,j),true);
  }

  avro_msgp( "there were " + std::to_string(nb_fixed) + " fixed elements with a crust of " + std::to_string(crust.size()) + " elements" );
  mpi::barrier();
}

template<typename type>
void
fix_conforming( Topology<type>& topology , const std::vector<VertexMetric>& metrics ) {

  MetricAttachment attachment( topology.points() , metrics );
  MetricField<type> field(topology,attachment);

  std::vector<index_t> edges;
  topology.get_edges(edges);

  std::vector< std::vector<index_t> > v2v( topology.points().nb() );
  for (index_t k = 0; k < edges.size()/2; k++) {
    v2v[ edges[2*k]   ].push_back( edges[2*k+1] );
    v2v[ edges[2*k+1] ].push_back( edges[2*k] );
  }

  // a pass using this feature is not intended to produce a metric conforming mesh
  // but to clean up the interfaces that may be rejected and are much further from metric-conforming
  real_t lim = 3.0;//sqrt(2.0);
  for (index_t k = 0; k < v2v.size(); k++) {
    bool ok = true;

    for (index_t j = 0; j < v2v[k].size(); j++) {

      real_t lk = field.length( topology.points() , k , v2v[k][j] );
      if (lk > lim || lk < 1/lim) ok = false;
    }

    if (ok)
      topology.points().set_fixed(k,true);
  }

  index_t nb_fixed = 0;
  for (index_t k = 0; k < topology.points().nb(); k++) {
    if (topology.points().fixed(k)) nb_fixed++;
  }
  //printf("there are %lu fixed points out of %lu\n",nb_fixed,topology.points().nb());
}

void
synchronize_metrics( const Points& points , std::vector<VertexMetric>& metrics0 ) {

  avro_assert( points.nb() == metrics0.size() );
  avro_assert( metrics0.size() > 0 );

  coord_t n = metrics0[0].n();

  // gather all the metrics and their associated global indices
  std::vector<index_t> indices;
  std::map<index_t,index_t> global2local;
  std::vector<VertexMetric> metrics;
  for (index_t k = 0; k < points.nb(); k++) {

    if (points.fixed(k)) {
      indices.push_back( points.global(k) );
      metrics.push_back( metrics0[k] );
      global2local.insert( {points.global(k) , k } );
    }
  }

  index_t rank = mpi::rank();
  index_t nb_rank = mpi::size();
  std::map<index_t,std::vector<VertexMetric>> global_metrics;
  if (rank == 0) {

    // append the metrics from the root
    for (index_t k = 0; k < indices.size(); k++) {
      global_metrics[indices[k]].push_back(metrics[k]);
    }

    for (index_t k = 1; k < nb_rank; k++) {

      // wait to receive our metrics from the processor
      std::vector<index_t> indices_k = mpi::receive<std::vector<index_t>>( k , TAG_METRIC_INDICES );
      std::vector<real_t> metric_data_k = mpi::receive<std::vector<real_t>>( k , TAG_METRIC_ENTRIES );

      // append the metrics from the processor
      index_t j = 0;
      for (index_t i = 0; i < indices_k.size(); i++) {

        index_t idx = indices_k[i];

        VertexMetric m(n);
        for (index_t r = 0; r < n; r++)
        for (index_t c = r; c < n; c++)
          m(r,c) = metric_data_k[j++];

        global_metrics[idx].push_back(m);
      }
    }
  }
  else {
      std::vector<real_t> metric_data;
    for (index_t k = 0; k < indices.size(); k++)
    for (index_t i = 0; i < n; i++)
    for (index_t j = i; j < n; j++)
      metric_data.push_back( metrics[k](i,j) );

    mpi::send( mpi::blocking{} , indices , 0 , TAG_METRIC_INDICES );
    mpi::send( mpi::blocking{} , metric_data , 0 , TAG_METRIC_ENTRIES );
  }
  mpi::barrier();

  if (rank == 0) {

    // average the metrics and then send them out to the processors
    std::map<index_t,std::vector<VertexMetric>>::const_iterator it;
    for (it = global_metrics.begin(); it != global_metrics.end(); ++it) {

      const index_t idx = it->first;
      const std::vector<VertexMetric>& metrics_k = it->second;
      std::vector<real_t> alpha( metrics_k.size() , 1./metrics_k.size() );
      try {
        VertexMetric m = numerics::interp( alpha , metrics_k );
        avro_assert( numerics::det(m) > 0 );
        global_metrics[idx] = {m};
      }
      catch (...) {
        print_inline(alpha);
        for (index_t i = 0; i < metrics_k.size(); i++)
          metrics_k[i].display();
      }
    }

    // send the new metric to each processor
    std::vector<index_t> idx;
      std::vector<real_t> metric_data;
    for (it = global_metrics.begin(); it != global_metrics.end(); ++it) {
      index_t global = it->first;
      avro_assert( it->second.size() == 1 ); // this has been reduced to 1 now
      const VertexMetric& m = it->second[0];
      for (index_t i = 0; i < n; i++)
      for (index_t j = i; j < n; j++)
        metric_data.push_back( m(i,j) );
      idx.push_back( global );

      if (global2local.find(global) != global2local.end())
        metrics0[global2local[global]] = m;
    }

    for (index_t k = 1; k < nb_rank; k++) {
      mpi::send( mpi::blocking{} , idx , k , TAG_METRIC_INDICES );
      mpi::send( mpi::blocking{} , metric_data , k , TAG_METRIC_ENTRIES );
    }
  }
  else {

    // receive the averaged metric from the root processor and assign it into the original metrics array
    std::vector<index_t> idx = mpi::receive< std::vector<index_t> >( 0 , TAG_METRIC_INDICES );
    std::vector<real_t> metric_data = mpi::receive< std::vector<real_t> >( 0 , TAG_METRIC_ENTRIES );

    index_t N = n*(n+1)/2;
    for (index_t k = 0; k < idx.size(); k++) {

      index_t global = idx[k];
      if (global2local.find(global) == global2local.end()) continue;

      index_t j = 0;
      for (index_t r = 0; r < n; r++)
      for (index_t c = r; c < n; c++)
        metrics0[global2local.at(global)](r,c) = metric_data[ N*k + j++ ];
    }
  }
  mpi::barrier();
}

template<typename type>
void
add_mantle( Topology<type>& topology , std::vector<index_t>& crust ) {

  topology.build_structures();

  std::set<index_t> pts;
  for (index_t k = 0; k < crust.size(); k++) {
    for (index_t j = 0; j < topology.nv(crust[k]); j++)
      pts.insert( topology(crust[k],j) );
  }

  std::vector<index_t> ball;
  for (std::set<index_t>::iterator it = pts.begin(); it != pts.end(); ++it) {
    ball.clear();
    topology.inverse().ball( *it , ball );
    for (index_t j = 0; j < ball.size(); j++)
      crust.push_back( ball[j] );
  }
  uniquify(crust);
  mpi::barrier();
}

template<typename type>
void
AdaptationManager<type>::adapt() {

  coord_t dim = topology_.points().dim();
  coord_t number = topology_.number();

  // clear any previous processor-processor pairs that may have been stored
  pairs_.clear();

  // setup some parameters
  bool limit_metric = params_["limit metric"];
  index_t max_passes = params_["max parallel passes"];
  index_t nb_part = mpi::size();
  index_t pass = 0;
  bool done = false;
  bool freeze_conforming = false;
  index_t force_partition_count = params_["force partition count"];
  index_t debug_level = params_["debug level"];
  bool write_mesh = params_["write mesh"];

  clock_t TIMER = clock();
  time_process_     = 0.0;
  time_synchronize_ = 0.0;
  time_adapt_       = 0.0;
  time_migrate_     = 0.0;
  time_partition_   = 0.0;
  time_exchange_    = 0.0;

  // loop over number of passes
  while (!done) {

    if (rank_ == 0) printf("\n*** pass %lu ***\n\n",pass);
    mpi::barrier();

    // fix the boundary of the topology
    // and move the fixed points to the beginning of the points structure
    if (rank_ == 0) TIMER = clock();
    if (rank_ == 0) printf("--> fixing partition boundaries\n");
    mpi::barrier();
    fix_boundary();
    mpi::barrier();

    // option to limit the metrics
    try {
      // only limit metrics if this is the first pass (if limiting is used at all)
      if (pass == 0 && params_["limit metric"]) {
        MetricAttachment attachment( topology_.points() , metrics_ );
        attachment.limit(topology_,2.0,true);
        for (index_t k = 0; k < attachment.nb(); k++)
          metrics_[k] = attachment[k];

        // we need to synchronize the metrics across processors that may have been modified by metric limiting
        synchronize_metrics( topology_.points() , metrics_ );
      }
      // do not limit the metrics in the call to the mesh adaptation below
      params_.set("limit metric", false);
    }
    catch(...) {
      printf("pass %lu: there was an error limiting the metric on processor %lu :(\n",pass,rank_);
      mpi::abort(1);
    }

    // check if we want to freeze points in the topology whose surrounding edge lengths are metric conforming
    // this gives a bit of a performance boost, but should be used with caution
    if (freeze_conforming) fix_conforming( topology_ , metrics_ );

    // create the mesh we will write to
    Mesh mesh(number,dim);
    std::shared_ptr<Topology<type>> ptopology = std::make_shared<Topology<type>>(mesh.points(),number);
    ptopology->TopologyBase::copy(topology_);
    mesh.add(ptopology);
    topology_.points().copy(mesh.points());
    Mesh mesh_out(number,dim);

    // option to write out the input mesh (for debugging)
    library::meshb writer;
    if (number < 4 && debug_level > 0)
      writer.write(mesh,"input-proc"+std::to_string(rank_)+".mesh",false);

    // setup the adaptation
    params_.set("output redirect" , "adaptation-output-proc"+std::to_string(rank_)+".txt" );
    params_.set("swapout" , false );
    params_.set("export boundary", false);
    params_.set("prefix" , "processor"+std::to_string(rank_)+"_pass"+std::to_string(pass) );
    params_.set("write mesh" , false ); // we don't want each processor to write the mesh it generates
    AdaptationProblem problem = {mesh,metrics_,params_,mesh_out};
    if (rank_ == 0) time_process_ += real_t( clock() - TIMER )/real_t(CLOCKS_PER_SEC);
    bool abort = false;
    try {
      // do the adaptation!
      clock_t t0 = 0, t1;
      if (rank_ == 0) {
        printf("--> adapting mesh on %lu partitions...\n",nb_part);
        t0 = clock();
      }
      ::avro::adapt<type>( problem );
      mpi::barrier();
      if (rank_ == 0) {
        t1 = clock();
        printf("--> time = %4.3g seconds.\n",real_t(t1-t0)/real_t(CLOCKS_PER_SEC));
        time_adapt_ += real_t(t1-t0)/real_t(CLOCKS_PER_SEC);
      }

      // clear the topology and copy in the output topology
      if (rank_ == 0) TIMER = clock();
      topology_.clear();
      topology_.points().clear();
      topology_.TopologyBase::copy( mesh_out.template retrieve<type>(0) );
      mesh_out.points().copy( topology_.points() );

      // unset any global indices in unfixed points, since these need to be renumbered
      for (index_t k = 0; k < topology_.points().nb(); k++) {
        if (topology_.points().fixed(k)) continue;
        topology_.points().set_global( k , 0 );
      }
      if (rank_ == 0) time_process_ += real_t(clock() - TIMER)/real_t(CLOCKS_PER_SEC);
    }
    catch(...) {
      printf("there was an error adapting the mesh on processor %lu :(\n",rank_);
      abort = true;
    }
    mpi::barrier();
    if (abort) {
	    mpi::abort(1);
    }

    // synchronize all the global point indices
    if (rank_ == 0) TIMER = clock();
    synchronize();
    if (rank_ == 0) time_synchronize_ += real_t( clock() - TIMER) / real_t(CLOCKS_PER_SEC);

    // option to retrieve and export the full mesh
    // we do not count this towards the time
    Mesh m(number,dim);
    std::shared_ptr<Topology<type>> topology_out = std::make_shared<Topology<type>>(m.points(),number);
    retrieve(*topology_out.get());
    if (rank_ == 0 && write_mesh) {

      index_t adapt_iter = params_["adapt iter"];

      nlohmann::json jm;
      jm["type"]     = topology_out->type_name();
      jm["dim"]      = topology_out->points().dim();
      jm["number"]   = topology_out->number();
      jm["vertices"] = topology_out->points().data();
      jm["elements"] = topology_out->data();
      jm["geometry"] = params_["geometry"];
      jm["nb_ghost"] = topology_out->points().nb_ghost();

      std::ofstream file("mesh-adapt"+std::to_string(adapt_iter)+"-pass"+std::to_string(pass)+".avro");
      file << jm;
    }

    #if 0
    // save the metric data for this processor
    // the metrics have been migrated so each processor own the metrics for the vertices it needs
    // we need to save the global id for vertex
    std::vector<real_t> metric_data( nb_metric_rank * metrics_.size() );
    index_t count = 0;
    for (index_t k = 0; k < metrics_.size(); k++)
    for (index_t j = 0; j < nb_metric_rank; j++)
     metric_data[count++] = metrics[k].data(j);

     avro_assert( metrics.size() == manager.topology().points().nb() );
     std::vector<index_t> global( metrics.size() );
     for (index_t k = 0; k < metrics.size(); k++)
       global[k] = manager.topology().points().global(k);

     nlohmann::json jm;
     jm["metrics"] = metric_data;
     jm["global"]  = global;

     std::ofstream file("metric-proc"+std::to_string(rank)+"-iter"+std::to_string(iter)+".metric");
     file << jm;
     file.close();

    #endif

    // report the min and max volumes in the mesh (this helps determine if something bad is happening...)
    if (rank_ == 0) {
      std::vector<real_t> volumes( topology_out->nb() );
      real_t v0 = topology_out->volume();
      topology_out->get_volumes(volumes);
      real_t min_v0 = *std::min_element( volumes.begin() , volumes.end() );
      real_t max_v0 = *std::max_element( volumes.begin() , volumes.end() );
      printf("--> volume = %g, min = %g, max = %g\n",v0,min_v0,max_v0);
    }

    // check if the user is forcing a number of partitions
    if (force_partition_count == 0) {

      // how many elements should we target per processor?
      index_t nb_elem_per_core = params_["elems per processor"];

      index_t nb_elem = topology_out->nb();
      if (rank_ == 0) {
        for (index_t k = 0; k < index_t(mpi::size()); k++)
          mpi::send( mpi::blocking{} , int(nb_elem) , k , TAG_MISC );
      }
      else {
        nb_elem = mpi::receive<int>( 0 , TAG_MISC );
      }

      nb_part = mpi::size();
      freeze_conforming = false;
      if (pass == max_passes-1 || max_passes <= 2) {
        done = true;
      }
      else if (pass == max_passes-2 && params_["allow serial"]) {
        freeze_conforming = true;
        nb_part = 1;
      }
      else {
        nb_part = index_t(nb_elem/nb_elem_per_core);
        if (nb_part == 0) nb_part = 1;
      }
      if (nb_part > index_t(mpi::size())) nb_part = mpi::size();
    }
    else {
      nb_part = force_partition_count;
      if (pass == max_passes-1) done = true;
    }

    // balance the mesh for the next pass, or for the caller
    if (rank_ == 0) TIMER = clock();
    if (rank_ == 0) printf("--> performing load balance & interface migration for %lu elements on %lu partitions\n",topology_out->nb(),nb_part);
    migrate_balance( nb_part );
    if (rank_ == 0) time_migrate_ += real_t( clock() - TIMER ) / real_t(CLOCKS_PER_SEC);
    pass++;
  }

  // reset some parameters to their original values
  params_.set("limit metric" , limit_metric );
  params_.set("write mesh" , write_mesh );
}

template<typename type>
class PartitionField : public Field<type,real_t> {
public:
  PartitionField( const Topology<type>& topology ) :
    Field<type,real_t>(topology,1,DISCONTINUOUS) {
    this->element().set_basis( BasisFunctionCategory_Lagrange );
  }

  index_t nb_rank() const override { return 1; }
  std::string get_name( index_t j ) const override {
    return std::to_string(j);
  }
};

template<typename type>
index_t
AdaptationManager<type>::get_total_nb_points() const {
  index_t nb_rank = mpi::size();

  // get the maximum global index
  index_t nb_points = 0;
  for (index_t k = 0; k < topology_.points().nb(); k++) {
    if (topology_.points().global(k) > nb_points)
      nb_points = topology_.points().global(k);
  }

  if (rank_ == 0) {
    for (index_t k = 1; k < nb_rank; k++) {
      index_t nb_points_k = index_t(mpi::receive<int>( k , TAG_MISC ));
      nb_points = std::max(nb_points,nb_points_k);
    }
  }
  else {
    mpi::send( mpi::blocking{} , int(nb_points) , 0 , TAG_MISC );
  }
  mpi::barrier();
  return nb_points;
}

template<typename type>
void
AdaptationManager<type>::append_partition( index_t p , Topology<type>& topology , Topology_Partition<type>& partition ,
                                           std::vector<bool>& created , std::vector<real_t>& partition_index ) const {
  for (index_t j = 0; j < partition.points().nb(); j++) {

    // retrieve the global index of this point
    index_t g = partition.points().global(j);
    if (created[g]) {
      // this point was already created
      continue;
    }
    created[g] = true;
    for (index_t d = 0; d < topology.points().dim(); d++) {
      topology.points()[g][d] = partition.points()[j][d];
    }
    topology.points().set_entity(g,partition.points().entity(j));
  }

  for (index_t j=0;j<partition.nb();j++) {
    partition_index.push_back(p);
    std::vector<index_t> s = partition.get(j);
    for (index_t i=0;i<s.size();i++) {
      s[i] = partition.points().global(s[i]);
      avro_assert( s[i] < topology.points().nb() );
    }
    topology.add(s.data(),s.size());
  }
}

template<typename type>
void
AdaptationManager<type>::retrieve( Topology<type>& topology ) {
  // fill the topology with the elements
  // fill the points with the vertices
  // assign the globals into the points
  index_t nb_rank = mpi::size();
  index_t nb_points = get_total_nb_points() + 1; // +1 because of 0-bias
  if (rank_ == 0) {
    std::vector<bool> created(nb_points,false);

    // initialize all the points (data filled below)
    std::vector<real_t> zero(topology.points().dim(),0.0);
    for (index_t j = 0; j < nb_points; j++)
      topology.points().create(zero.data());

    // append the root partition into the output topology
    std::vector<real_t> partition_index;
    append_partition( 0 , topology , topology_ , created , partition_index );

    for (index_t k = 1; k < nb_rank; k++) {

      // receive the partition from processor k and append it into the output topology
      Topology_Partition<type> partition( topology_.points().dim() , topology_.points().udim() , topology_.number() );
      partition.set_entities( topology_.entities() );
      partition.receive(k);
      append_partition( k , topology , partition , created , partition_index );
    }

    // fill the partition field
    if (topology.nb() < 100000) {
      // field construction will be super slow if there are too many elements
      field_ = std::make_shared<PartitionField<type>>(topology);
      topology.element().set_basis( BasisFunctionCategory_Lagrange );
      field_->build();
      for (index_t k = 0; k < topology.nb(); k++)
      for (index_t j = 0; j < topology.nv(k); j++)
        (*field_)(k,j) = partition_index[k];
      topology.fields().make("partition",field_);
    }
  }
  else {
    // send the working topology to the root processor
    topology_.send(0);
  }
  mpi::barrier();
}

template<typename type>
void
AdaptationManager<type>::reassign_metrics( const std::vector<VertexMetric>& metrics ) {
  metrics_ = metrics;
}

template class AdaptationManager<Simplex>;

#endif

} // avro
