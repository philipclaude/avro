//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

namespace avro
{

template<>
Topology<Polytope>::Topology( Points& points , coord_t number , coord_t order ) :
  TopologyBase(points,number,TableLayout_Jagged,"polytope"),
  element_( number , order , points.incidence() ),
  neighbours_(*this),
  inverse_(*this)
{}

template<>
void
Topology<Polytope>::orient( index_t* v , const index_t nv , real_t* p )
{
  avro_implement;
}

template class Topology<Polytope>;
template class Tree<Topology<Polytope>>;
template void Topology<Polytope>::construct( std::shared_ptr<Topology<Simplex>>& node , Topology<Simplex>& ) const;
template void Topology<Polytope>::construct( std::shared_ptr<Topology<Polytope>>& node , Topology<Polytope>& ) const;

} // avro
