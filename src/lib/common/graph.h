//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_COMMON_GRAPH_H_
#define avro_LIB_COMMON_GRAPH_H_

namespace avro
{

template<typename Derived_t>
class Graph
{
public:

  void path( const Node_t& node_i , const Node_t& node_j ) const;

  real_t distance( const Node_t& node_i , const Node_t& node_j ) const;

private:
  Derived_t& derived() { return *static_cast<Derived_t*>(this); }
  const Derived& derived() { return *static_cast<const Derived_t*>(this); }

};

} // avro

#endif
