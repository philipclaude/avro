#ifndef avro_LIB_ADAPTATION_CAVITY_H_
#define avro_LIB_ADAPTATION_CAVITY_H_

#include "common/error.h"
#include "common/types.h"

#include "adaptation/implied_metric.h"

#include <map>
#include <string>
#include <vector>

namespace avro
{

class Entity;

void
geometryParams( Entity* e , const Points& points , const index_t* v , const index_t nv , real_t* params );

template<typename type>
class Cavity : public Topology<type>
{

public:
  Cavity( Topology<type>& _topology );

  // computation of the cavity given a point index and associated coordinate
  // initial cavity provided in C0, initialized be specialized operators
  bool compute( const index_t p , real_t* x , const std::vector<index_t>& C0 );
  bool enlarge( bool verbose=false );
  bool findGeometry( real_t* x , std::vector<index_t>& C0 );
  void apply();

  // computation of the nodes in the cavity and whether any nodes are
  // removed by the application of the operator
  // TODO refactor naming to "vertices" not "nodes" to be consistent
  // with my terminology (not that of coupez)
  void computeNodes();
  void computeRemovedNodes();
  index_t nb_removedNodes() const { return removedNodes_.size(); }
  const std::vector<index_t>& nodes() const { return nodes_; }
  index_t nb_cavity() const { return cavity_.size(); }
  const std::vector<index_t>& cavity() const { return cavity_; }
  index_t nb_bnd() const { return boundary_.nb(); }
  index_t nb_insert() const { return Topology<type>::nb(); }
  index_t nb_nodes() const { return nodes_.size(); }

  bool removesNodes() const { return nb_removedNodes()>0; }
  index_t removedNode( const index_t k ) const { return removedNodes_[k]; }
  index_t cavity( const index_t k ) const { return cavity_[k]; }
  bool contains( const index_t c ) const;

  // computation of the cavity boundary and the geometry topology
  bool computeBoundary();
  Topology<type>& boundary() { return boundary_; }
  void computeGeometry( Entity* g , Topology<type>& G , std::map<index_t,index_t>& v2u , std::vector<index_t>& u2v );

  real_t minvol() const { return minvol_; }

  bool& rethrow() { return rethrow_; }

  bool positiveImpliedMetrics();

  bool philipcondition();

  void print() const;

  void add( const index_t* v , const index_t nv );
  void addNode( const index_t node );
  void addRemovedNode( const index_t node );
  void addCavity( const index_t elem );

  void clear();

  // used to store the insertion values
  std::vector<index_t>& inserted() { return inserted_; }
  std::vector<bool>& removes() { return removes_; }

  void setEnlarge( bool x ) { enlarge_ = x; }
  void setIgnore( bool x ) { ignore_ = x; }

  const Topology<type>& topology() const { return topology_; }

  real_t& sign() { return sign_; }

protected:
  Topology<type>& topology_;
  bool nodeRemovalAllowed_;
  bool enlarge_;
  bool checkVisibility_;
  std::string info_;

private:
  index_t star_;
  std::vector<real_t> point_;

  Topology<type> boundary_; // boundary of the elements
  std::vector<index_t> cavity_;  // elements
  std::vector<index_t> nodes_; // nodes of the cavity elements
  std::vector<index_t> removedNodes_; // deleted nodes
  std::vector<index_t> idx_; // index of the star for each boundary facet
  std::vector<index_t> inserted_; // index of the insertion in the topology
  std::vector<bool> removes_; // whether the cavity element is actually removed from the topology

  // minimum volume the cavity elements can create
  real_t minvol_;
  real_t sign_;

  // some error-tracking/debugging stuff
  index_t nb_error_;
  bool rethrow_;

  // element-implied metric calculator
  ElementImpliedMetric<type> mk_;

  bool ignore_;

};

} // avro

#endif
