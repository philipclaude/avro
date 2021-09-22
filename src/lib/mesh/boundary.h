//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_MESH_BOUNDARY_H_
#define avro_MESH_BOUNDARY_H_

#include "avro_types.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <map>

namespace avro
{

namespace BoundaryUtils
{
  Entity*
  geometryFacet( const Points& points , const index_t* v , index_t nv , bool elem=false );
}

class Entity;

class BoundaryPoints : public Points
{
public:
  BoundaryPoints( const Points& points , bool meta = true );

  index_t index( index_t k ) const { return local2global_.at(k); }

private:
  std::map<index_t,index_t> local2global_;
};

template<typename type>
class Boundary : public Topology<type>
{
public:
  Boundary( const Topology<type>& _topology ) :
      Topology<type>( _topology.points() , _topology.number() ),
      topology_(_topology),
      nodes_(topology_.points(),0),
      edges_(topology_.points(),1),
      triangles_(topology_.points(),2),
      tetrahedra_(topology_.points(),3)
  {}

  void extractall();
  void extract(bool interior=false);

  real_t volume( const index_t k )
  {
    real_t vol = 0.;
    Topology<type>& t = this->child(k);
    for (index_t j=0;j<t.nb();j++)
    {
      std::vector<const real_t*> xk;
      t.get_elem(j,xk);
      vol += numerics::volume( xk , t.points().dim() );
    }
    return vol;
  }

  void add( index_t* v , index_t nv , Entity* e );

  index_t indexof( Entity* e );

  // returns the entity attached to the k'th boundary topology
  Entity* entity( const index_t k );
  const std::vector<Entity*>& entities() const { return entity_; }

  bool check() const;

  void print() const;

  const Topology<Simplex> nodes() const { return nodes_; }
  const Topology<Simplex> edges() const { return edges_; }
  const Topology<Simplex> triangles() const { return triangles_; }
  const Topology<Simplex> tetrahedra() const { return tetrahedra_; }

  const std::vector<index_t>& node_entities() const { return node_entities_; }
  const std::vector<index_t>& edge_entities() const { return edge_entities_; }
  const std::vector<index_t>& triangle_entities() const { return triangle_entities_; }
  const std::vector<index_t>& tetrahedron_entities() const { return tetrahedron_entities_; }

private:
  const Topology<type>& topology_;

  std::map< Entity* , index_t > entity2child_;
  std::map< index_t , index_t > bnd2child_;

  std::vector<Entity*> entity_;

  Topology<type> nodes_;
  Topology<type> edges_;
  Topology<type> triangles_;
  Topology<type> tetrahedra_;

  std::vector<index_t> node_entities_;
  std::vector<index_t> edge_entities_;
  std::vector<index_t> triangle_entities_;
  std::vector<index_t> tetrahedron_entities_;

};


} // avro

#endif
