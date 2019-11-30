#ifndef LUNA_LIB_MESH_H_
#define LUNA_LIB_MESH_H_

#include "common/types.h"

#include "mesh/topology.h"

#include <memory>
#include <vector>

namespace luna
{

template<typename type>
class Mesh
{
public:
  typedef Topology<type> Topology_t;
  typedef std::shared_ptr<Topology_t> Topology_ptr;

  index_t nb_topologies() const { return topology_.size(); }

  void add( Topology_ptr topology );
  
  Topology_t& topology( index_t k) { return *topology_[k].get(); }
  const Topology_t& topology( index_t k ) const { return *topology_[k].get(); }

private:
  std::vector<Topology_ptr> topology_;
};

} // luna

#endif
