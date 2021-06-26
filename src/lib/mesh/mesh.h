//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MESH_H_
#define avro_LIB_MESH_H_

#include "types.h"

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
  void set_number( coord_t number ) { number_ = number; }

  void add( Topology_ptr topology ) { topology_.push_back(topology); }

  TopologyBase& topology(index_t k) { return *topology_[k].get(); }
  const TopologyBase& topology(index_t k) const { return *topology_[k].get(); }

  Topology_ptr topology_ptr( index_t k ) { return topology_[k]; }

  template<typename type>
  Topology<type>&
  retrieve( index_t k )
  {
    avro_assert_msg( k < nb_topologies() , "requested topology %lu of %lu" , k , nb_topologies() );
    avro_assert_msg( topology_[k]->type_name() == type::type_name()  ,
      "requested topology of type %s but have %s" , topology_[k]->type_name().c_str() , type::type_name().c_str() );
    return static_cast<Topology<type>&>(*topology_[k].get());
  }

  template<typename type>
  const Topology<type>&
  retrieve( index_t k ) const
  {
    avro_assert_msg( k < nb_topologies() , "requested topology %lu of %lu" , k , nb_topologies() );
    avro_assert_msg( topology_[k]->type_name() == type::type_name()  ,
      "requested topology of type %s but have %s" , topology_[k]->type_name().c_str() , type::type_name().c_str() );
    return static_cast<Topology<type>&>(*topology_[k].get());
  }

  template<typename type>
  std::shared_ptr<Topology<type>>
  retrieve_ptr( index_t k )
  {
    avro_assert_msg( topology_[k]->type_name() == type::type_name()  ,
      "requested topology of type %s but have %s" , topology_[k]->type_name().c_str() , type::type_name().c_str() );
    return std::static_pointer_cast<Topology<type>>(topology_[k]);
  }

  Points& points() { return points_; }
  const Points& points() const { return points_; }

  void retrieve( std::vector<const TopologyBase*>& topologies ) const
  {
    topologies.clear();
    for (index_t k=0;k<topology_.size();k++)
    {
      topologies.push_back( static_cast<const TopologyBase*>(topology_[k].get()) );
      if (topology_[k]->type_name()=="simplex")
      {
        static_cast<const Topology<Simplex>&>(*topology_[k]).get_children_typed<TopologyBase>(topologies);
      }
      else
        avro_implement;
    }
  }

  template<typename type>
  void retrieve( std::vector<Topology<type>*>& topologies )
  {
    topologies.clear();
    for (index_t k=0;k<topology_.size();k++)
    {
      if (topology_[k]->type_name()==type::type_name())
      {
        topologies.push_back( static_cast<Topology<type>*>(topology_[k].get()) );
        static_cast<Topology<type>&>(*topology_[k]).get_children(topologies);
      }
      else
        avro_implement;
    }
  }

  void retrieve( std::vector<const TopologyBase*>& topologies )
  {
    topologies.clear();
    for (index_t k=0;k<topology_.size();k++)
    {
      topologies.push_back( topology_[k].get() );
      topology_[k]->get_topologies(topologies);
    }
  }

  template<typename type>
  void retrieve( Topology<type>& topology )
  {
    for (index_t k=0;k<topology_.size();k++)
    {
      if (topology_[k]->type_name()==type::type_name())
        static_cast<Topology<type>&>(*topology_[k]).get_elements(topology);
    }
  }

protected:
  Points points_;
  coord_t number_;
  std::vector<Topology_ptr> topology_;
};



} // avro

#endif
