#ifndef avro_LIB_MESH_H_
#define avro_LIB_MESH_H_

#include "common/types.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include <memory>
#include <vector>

namespace avro
{

class Mesh
{
protected:
  typedef std::shared_ptr<TopologyBase> Topology_ptr;

public:
  Mesh( coord_t dim );
  Mesh( coord_t number , coord_t dim );

  index_t nb_topologies() const { return topology_.size(); }

  coord_t number() const { return number_; }

  void add( Topology_ptr topology ) { topology_.push_back(topology); }

  TopologyBase& topology(index_t k) { return *topology_[k].get(); }
  const TopologyBase& topology(index_t k) const { return *topology_[k].get(); }

  template<typename type>
  Topology<type>&
  retrieve( index_t k )
  {
    avro_assert_msg( topology_[k]->type_name() == type::type_name()  ,
      "requested topology of type %s but have %s" , topology_[k]->type_name().c_str() , type::type_name().c_str() );
    return static_cast<Topology<type>&>(*topology_[k].get());
  }

  template<typename type>
  const Topology<type>&
  retrieve( index_t k ) const
  {
    avro_assert_msg( topology_[k]->type_name() == type::type_name()  ,
      "requested topology of type %s but have %s" , topology_[k]->type_name().c_str() , type::type_name().c_str() );
    return static_cast<Topology<type>&>(*topology_[k].get());
  }

  Points& points() { return points_; }
  const Points& points() const { return points_; }

protected:
  Points points_;
  coord_t number_;
  std::vector<Topology_ptr> topology_;
};

} // avro

#endif
