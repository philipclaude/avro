//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MASTER_POLYTOPE_H_
#define avro_LIB_MASTER_POLYTOPE_H_

#include "common/table.h"

#include "element/element.h"
#include "element/simplex.h"

#include <set>

namespace avro
{

template<typename type> class Topology;
template<typename type> class SimplicialDecomposition;

class Polytope : public Element<Polytope>
{

public:
  Polytope( coord_t number , coord_t order , const Table<int>& incidence );
  Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence );

  static std::string type_name() { return "polytope"; }
  static TableLayoutCategory layout() { return TableLayout_Jagged; }

  const Table<int>& incidence() const { return incidence_; }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

  real_t volume( const Points& points , const index_t* v , index_t nv ) const;
  real_t volume( const std::vector<const real_t*>& x , const coord_t dim ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const
    { avro_implement; }

    // for switching between the hrep and vrep
  void hrep( const index_t* v , index_t nv , std::vector<int>& facets ) const;
  void vrep( const index_t* v , index_t nv , int facet , std::vector<index_t>& points ) const;

  // edge retrieval
  void edges( const index_t* v , index_t nv , std::vector<index_t>& e  ) const;
  bool is_edge( index_t v0 , index_t v1 ) const;
  bool is_edge( const std::vector<index_t>& b0 , const std::vector<index_t>& b1 ) const;
  bool is_edge( const std::vector<int>& b0 , const std::vector<int>& b1 ) const;
  index_t edge( index_t k , index_t j ) const
    { avro_assert_not_reached; return 0; }

  std::vector<index_t> triangulate( const index_t* v , index_t nv , SimplicialDecomposition<Polytope>& triangulation , index_t parent , std::set<int>& h ) const;

  bool& fullmesh() { return fullmesh_; }

private:
  void triangulate( coord_t number , Table<index_t>& table , Points& points , const index_t* v , index_t nv ) const;

  Simplex simplex_;
  const Table<int>& incidence_;

  bool fullmesh_; // whether the  polytope is responsible for a full mesh
};

} // avro

#endif
