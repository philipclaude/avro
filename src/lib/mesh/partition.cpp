#include "common/process.h"
#include "common/tools.h"

#include "element/simplex.h"

#include "mesh/partition.h"
#include "mesh/topology.h"

#include <set>

#ifdef AVRO_MPI

#include "common/mpi.hpp"

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
Partition<type>::get( std::vector< std::shared_ptr<Topology<type>>>& parts ) const
{
}

template class Partition<Simplex>;

} // avro

#endif
