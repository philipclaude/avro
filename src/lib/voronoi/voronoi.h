//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_MESH_DELAUNAY_VORONOI_H_
#define AVRO_MESH_DELAUNAY_VORONOI_H_

#include "common/nearest_neighbours.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/field.hpp"
#include "mesh/mesh.h"
#include "mesh/topology.h"

#include "numerics/predicates.h"

namespace avro
{

class Delaunay;
class RestrictedVoronoiSimplex;
class Entity;

namespace delaunay
{

typedef struct
{
  std::vector<int> indices;
} SymbolicVertex;

// needed to create a set/map of elements
inline bool
operator==( const SymbolicVertex& fx , const SymbolicVertex& fy )
{
  // assumes fx and fy have the same topological dimension
  // and that the indices are sorted
  avro_assert( fx.indices.size()==fy.indices.size() );
  for (index_t j=0;j<fx.indices.size();j++)
    if (fx.indices[j]!=fy.indices[j])
      return false;
  return true;
}

// needed to create a map of elements
inline bool
operator<( const SymbolicVertex& f , const SymbolicVertex& g )
{
  // lexicographically compare the indices
  return std::lexicographical_compare(f.indices.begin(), f.indices.end(),
                                      g.indices.begin(), g.indices.end());
}

class RVDFacets
{
public:
  RVDFacets( const Topology<Simplex>& topology );

  void create();
  int facet( const std::vector<index_t>& f ) const;
  void print() const;

private:
  std::string generate( const std::vector<index_t>& f ) const;
  int lookup( const std::string& s , int& id ) const;

  std::map<std::string,int> store_;
  const Topology<Simplex>& topology_;

};

class VoronoiSites : public Field<Polytope,real_t>
{
public:
  VoronoiSites( Topology<Polytope>& rvd ) :
    Field<Polytope,real_t>(rvd,0,DISCONTINUOUS)
  {}

  void add_cell( const real_t z )
  {
    //Field<real_t>::add(z+1);
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("sites"); return result;}
};

} // delaunay

} // avro

#endif
