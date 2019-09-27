// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_MESH_TOPOLOGY_H_
#define URSA_MESH_TOPOLOGY_H_

#include "common/data.h"
#include "common/json.h"
#include "common/tree.h"
#include "common/types.h"

#include "master/master.h"

#include <vector>

namespace ursa
{

class Vertices;
class ClippingPlane;

class TopologyBase : public Data<index_t>, public TreeNodeBase
{
public:
  TopologyBase( Vertices& vertices , const coord_t number ) :
    Data<index_t>(true) , vertices_(vertices) , number_(number)
  {}

  virtual ~TopologyBase() {}

  Vertices& vertices() const { return vertices_; }

  coord_t number() const { return number_; }
  void set_number( const coord_t _number ) { number_ = _number; }

  // virtual functions for leaf
  void copy( TopologyBase& topology );

  // index/cell retrieval
  std::vector<index_t> get( const index_t k ) const
    { return Data::get(k); }
  index_t operator() ( const index_t k , const index_t j ) const
    { return Data<index_t>::operator()(k,j); }
  index_t& operator() ( const index_t k , const index_t j )
    { return Data<index_t>::operator()(k,j); }

  virtual void getChildren( std::vector<TopologyBase*>& children ) const = 0;
  virtual void getPoints( std::vector<index_t>& p ) const = 0;
  virtual void getEdges( std::vector<index_t>& e , ClippingPlane* plane ) const = 0;
  virtual void getTriangles( std::vector<index_t>& t , Vertices& v, std::vector<index_t>& parent , ClippingPlane* plane ) const = 0;

  const std::string& name() const { return name_; }
  void setName( const std::string& _name ) { name_ = _name; }

  void setDummy( bool x ) { dummy_ = x; }
  bool dummy() const { return dummy_; }

protected:
  Vertices& vertices_;
  coord_t number_;
  bool dummy_;
  std::string name_;
};


template <typename type>
class Topology : public Tree< Topology<type> >, public TopologyBase
{

public:
  typedef std::shared_ptr<Topology<type>> Topology_ptr;

  Topology( Vertices& _vertices , const coord_t _number );
  Topology( Vertices& _vertices , const coord_t _number , const coord_t _order );
  Topology( Vertices& _vertices , const json& J );
};


template<typename Basis>
class Topology< Master<Simplex<Basis>> >
{
  using Tree< Topology<Master<Simplex<Basis>>> >::nb_children;
};

}

#endif
