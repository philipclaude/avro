#ifndef AVRO_LIB_MESH_PARTITION_H_
#define AVRO_LIB_MESH_PARTITION_H_

#include "common/types.h"

#include "mesh/points.h"

#include <map>
#include <memory>
#include <vector>

namespace avro
{

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

  void move_to_front( const std::vector<index_t>& pts );

  index_t local2global( index_t k ) const;
  index_t global2local( index_t k ) const;

private:
  Points points_;
  std::vector<index_t>      local2global_;
  std::map<index_t,index_t> global2local_;
  std::vector<Entity*> entities_;
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
