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
  void send( mpi::communicator& comm , index_t receiver ) const;
  void receive( mpi::communicator& comm , index_t sender );
  #endif

  void extract_points( const Points& points );
  void convert();

  index_t local2global( index_t k ) const;
  index_t global2local( index_t k ) const;

  void reset_indices();
  void set_local2global( index_t k , index_t global );
  void set_global2local( index_t global , index_t k );

  void map_indices( const std::map<index_t,index_t>& idx );
  void remove_indices( const std::vector<index_t>& idx );

  void compute_crust(); // should be removed
  void compute_crust( const std::vector<index_t>& halo , std::vector<index_t>& crust ) const;
  void compute_mantle( std::vector<index_t>& interior , std::vector<index_t>& exterior , std::vector<index_t>& mantle , const std::set<index_t>& halo ) const;

  const std::vector<index_t>& mantle() const { return mantle_; }
  const std::vector<index_t>& crust() const { return crust_; }
  const std::vector<index_t>& halo() const { return halo_; }

  bool check( const Points& points ) const;

protected:
  Points points_;
  std::vector<index_t>      local2global_;
  std::map<index_t,index_t> global2local_;
  std::vector<Entity*> entities_;
  std::vector<index_t> halo_;
  std::vector<index_t> crust_;
  std::vector<index_t> mantle_;
};

template<typename type>
class MergedPartitions : public Topology_Partition<type>
{
public:
  #ifdef AVRO_MPI
  void receive_partition( mpi::communicator& comm , index_t receiver );
  void receive_halo( mpi::communicator& comm , index_t receiver );
  #endif

private:
  std::vector<index_t> boundary_;
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
