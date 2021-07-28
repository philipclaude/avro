//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_CAVITY_H_
#define avro_LIB_ADAPTATION_CAVITY_H_

#include "common/error.h"
#include "avro_types.h"

#include "adaptation/implied_metric.h"

#include <map>
#include <string>
#include <unordered_set>
#include <vector>

namespace avro
{

class Entity;

void
geometry_params( Entity* e , const Points& points , const index_t* v , const index_t nv , real_t* params );

template<typename type>
class Cavity : public Topology<type>
{

public:
  Cavity( Topology<type>& _topology );

  // computation of the cavity given a point index and associated coordinate
  // initial cavity provided in C0, initialized be specialized operators
  bool compute( const index_t p , real_t* x , const std::vector<index_t>& C0 );
  bool enlarge( bool verbose=false );
  bool find_geometry( real_t* x , std::vector<index_t>& C0 );
  void apply();

  // computation of the nodes in the cavity and whether any nodes are
  // removed by the application of the operator
  // TODO refactor naming to "vertices" not "nodes" to be consistent
  // with my terminology (not that of coupez)
  void compute_nodes();
  void compute_removed_nodes();
  index_t nb_removed_nodes() const { return removed_nodes_.size(); }
  const std::vector<index_t>& nodes() const { return nodes_; }
  index_t nb_cavity() const { return cavity_.size(); }
  const std::vector<index_t>& cavity() const { return cavity_; }
  index_t nb_bnd() const { return boundary_.nb(); }
  index_t nb_insert() const { return Topology<type>::nb(); }
  index_t nb_nodes() const { return nodes_.size(); }

  bool removes_nodes() const { return nb_removed_nodes()>0; }
  index_t removed_node( const index_t k ) const { return removed_nodes_[k]; }
  index_t cavity( const index_t k ) const { return cavity_[k]; }
  bool contains( const index_t c ) const;

  void copy( const Cavity<type>& cavity );

  // computation of the cavity boundary and the geometry topology
  bool compute_boundary();
  Topology<type>& boundary() { return boundary_; }
  void compute_geometry( Entity* g , Topology<type>& G , std::map<index_t,index_t>& v2u , std::vector<index_t>& u2v );

  real_t minvol() const { return minvol_; }

  bool& rethrow() { return rethrow_; }

  bool positive_implied_metrics();

  bool has_unique_elems();

  void print() const;

  void add( const index_t* v , const index_t nv );
  void add_node( const index_t node );
  void add_removed_node( const index_t node );
  void add_cavity( const index_t elem );

  void clear();

  // used to store the insertion values
  std::vector<index_t>& inserted() { return inserted_; }
  std::vector<bool>& removes() { return removes_; }

  void set_enlarge( bool x ) { enlarge_ = x; }
  void set_ignore( bool x ) { ignore_ = x; }
  void check_visibility( bool x ) { check_visibility_ = x; }

  const Topology<type>& topology() const { return topology_; }

  real_t& sign() { return sign_; }

  void set_entity( Entity* entity ) { entity_ = entity; }
  Entity* entity() { return entity_; }

  bool fixed() const;
  bool closed_boundary();

protected:
  Topology<type>& topology_;
  bool node_removal_allowed_;
  bool enlarge_;
  bool check_visibility_;
  std::string info_;

private:
  index_t star_;
  std::vector<real_t> point_;

  Topology<type> boundary_; // boundary of the elements
  std::vector<index_t> cavity_;  // elements
  std::vector<index_t> nodes_; // nodes of the cavity elements
  std::vector<index_t> removed_nodes_; // deleted nodes
  std::vector<index_t> idx_; // index of the star for each boundary facet
  std::vector<index_t> inserted_; // index of the insertion in the topology
  std::vector<bool> removes_; // whether the cavity element is actually removed from the topology

  std::unordered_set<index_t> cavity_set_;

  // minimum volume the cavity elements can create
  real_t minvol_;
  real_t sign_;

  // some error-tracking/debugging stuff
  index_t nb_error_;
  bool rethrow_;

  // element-implied metric calculator
  ElementImpliedMetric<type> mk_;

  bool ignore_;

protected:
  Entity* entity_;

};

} // avro

#endif
