#ifndef AVRO_LIB_MESH_PARTITION_H_
#define AVRO_LIB_MESH_PARTITION_H_

#include "common/types.h"

#include "mesh/points.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace avro
{

inline bool
fixed_facet( const std::vector<index_t>& facet , const Points& points )
{
  for (index_t j=0;j<facet.size();j++)
  {
    if (!points.fixed(facet[j]))
      return false;
  }
  return true;
}

template<typename type> class Topology;
class Entity;

// a partitioned topology that holds its own set of points
template<typename type>
class Topology_Partition : public Topology<type>
{
public:
  Topology_Partition( coord_t dim , coord_t udim , coord_t number );

  // set the entities we will use to look up geometry identifiers
  void set_entities( const std::vector<Entity*>& entities );
  Entity* lookup( int identifier , int number ) const;

  #ifdef AVRO_MPI
  void send( index_t receiver ) const;
  void receive( index_t sender );
  #endif

  void convert( const Points& points );

  bool check( const Points& points ) const;

protected:
  Points points_;
  std::vector<Entity*> entities_;
};

template<typename type>
class PartitionBoundary : public Topology<type>
{
public:
  PartitionBoundary( coord_t number );

  void send( index_t receiver ) const;
  void receive( index_t sender );

  void compute( const Topology<type>& topology , index_t partition );

  void add_left( const ElementIndices& f , index_t elemL , index_t partL , index_t elemL_global );
  void add_rite( const std::vector<index_t>& f , index_t elemL , index_t partL , index_t elemL_global );

  index_t elemL( index_t k ) const { return elemL_[k]; }
  index_t elemL_global( index_t k ) const { return elemL_global_[k]; }
  index_t partL( index_t k ) const { return partL_[k]; }

  index_t elemR( index_t k ) const { return elemR_[k]; }
  index_t elemR_global( index_t k ) const { return elemR_global_[k]; }
  index_t partR( index_t k ) const { return partR_[k]; }

  void ball( index_t p , std::vector<index_t>& B ) const;

  void append( const PartitionBoundary<type>& boundary );
  void fill( const PartitionBoundary<type>& interface );

  int find( const ElementIndices& f ) const;

  index_t nb() const { return facet_.size(); }

  const std::map<ElementIndices,index_t>& facets() const { return facet_; }

  void print() const;

private:
  Points dummy_; // should not be used, but we need this to inherit from topology so we can use the inverse
  std::vector<index_t> elemL_;
  std::vector<index_t> elemL_global_;
  std::vector<index_t> partL_;

  std::vector<index_t> elemR_;
  std::vector<index_t> elemR_global_;
  std::vector<index_t> partR_;

  std::map<ElementIndices,index_t> facet_;
};


template<typename type>
class Partition
{
public:
  Partition( const Topology<type>& topology );

  void initialize();
  void compute( index_t np );
  void get( std::vector< std::shared_ptr<Topology_Partition<type>> >& partitions ) const;

  const std::vector<index_t>& partition() const { return partition_; }

  void add_adjacency( index_t u , index_t v );

  bool weighted() const { return adjwgt_.size()==adjncy_.size(); }

  void compute_interface_points( std::vector<std::set<index_t>>& pts ) const;

  void print() const;

private:
  const Topology<type>& topology_;

  // coo storage format
  std::vector<real_t> weight_;
  std::vector<index_t> edges_;

  // csr storage format
  std::vector<index_t> xadj_;
  std::vector<index_t> adjncy_;
  std::vector<real_t> vwgt_; // vertex weights, size n
  std::vector<real_t> adjwgt_; // edge weights, size 2m

  std::vector<index_t> partition_;
};



} // avro

#endif
