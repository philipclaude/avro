//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_SPACETIME_H_
#define avro_LIB_LIBRARY_SPACETIME_H_

#include "mesh/mesh.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include <map>

namespace avro
{

class Entity;

template<typename type>
class Topology_Spacetime : public Topology<type>
{
public:
  Topology_Spacetime( const Topology<type>& topology );

  void extract();

  void write( const std::string& filename );
  void read( const std::string& filename );

  void clear();

  Entity* entity( Topology<type>* t ) const { return topology2entity_.at(t); }

private:
  Points points_;
  const Topology<type>& topology_;
  //Topology<type> slice_;
  std::map<Entity*,Topology<type>*> entity2topology_;
  std::map<Topology<type>*,Entity*> topology2entity_;
};

} // avro

#endif
