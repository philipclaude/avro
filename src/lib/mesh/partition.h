#ifndef AVRO_LIB_MESH_PARTITION_H_
#define AVRO_LIB_MESH_PARTITION_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace avro
{

template<typename type> class Topology;

template<typename type>
class Partition
{
public:
  Partition( const Topology<type>& topology );

  void initialize();
  void compute( index_t np );
  void get( std::vector< std::shared_ptr<Topology<type>> >& partitions ) const;

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
