#ifndef LUNA_MESH_TOPOLOGY_H_
#define LUNA_MESH_TOPOLOGY_H_

#include "common/table.h"
#include "common/json.h"
#include "common/tree.h"
#include "common/types.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/field.h"
#include "mesh/inverse.h"
#include "mesh/neighbours.h"

#include <string>
#include <vector>

namespace luna
{

class Points;
template<typename type> class Cavity;

class TopologyBase : public Table<index_t>
{
protected:
  TopologyBase( Points& vertices , const coord_t number , TableLayoutCategory category , const std::string& type ) :
    Table<index_t>(category,number+1),
    points_(vertices),
    number_(number),
    fields_(*this),
    type_(type)
  {}

public:
  virtual ~TopologyBase() {}

  Points& points() const { return points_; }

  coord_t number() const { return number_; }
  coord_t order() const { return order_; }

  void set_number( coord_t number ) { number_ = number; set_rank(number_+1); }

  void copy( TopologyBase& topology );

  virtual void get_points( std::vector<index_t>& p ) const = 0;
  virtual void get_edges( std::vector<index_t>& e ) const = 0;
  virtual void get_triangles( std::vector<index_t>& t ) const = 0;

  const std::string& name() const { return name_; }
  void setName( const std::string& _name ) { name_ = _name; }

  const std::string& type_name() const { return type_; }

  void set_dummy( bool x ) { dummy_ = x; }
  bool dummy() const { return dummy_; }

  template<typename derived_t>
  void attach( const std::string& name , std::shared_ptr<derived_t>& fld )
  {
    fields_.make( name , fld );
  }

  Fields& fields() { return fields_; }
  const Fields& fields() const { return fields_; }

protected:
  Points& points_;
  coord_t number_;
  coord_t order_;
  bool dummy_;
  std::string name_;
  Fields fields_;
  std::string type_;
};

template<typename type>
class Topology : public Tree<Topology<type>>, public TopologyBase
{
public:
  typedef Topology<type> Topology_t;
  typedef std::shared_ptr<Topology_t> Topology_ptr;
  using Tree<Topology_t>::nb_children;

  Topology( Points& vertices , coord_t number );
  Topology( Points& vertices , coord_t number , coord_t order );
  Topology( Points& vertices , const Topology<type>& topology , coord_t order );
  Topology( Points& vertices , const json& J );

  Topology_t& topology( index_t k ) { return Tree<Topology_t>::child(k); }

  type& master() { return master_; }
  const type& master() const { return master_; }

  coord_t number() const { return master_.number(); }
  coord_t order() const { return master_.order(); }

  bool ghost( index_t k ) const { luna_implement; return false; }
  index_t nb_real() const;
  index_t nb_ghost() const;

  void get_points( std::vector<index_t>& p ) const {}
  void get_edges( std::vector<index_t>& e ) const;
  void get_triangles( std::vector<index_t>& t ) const;

  void get_elem( index_t k , std::vector<real_t*>& X ) const;
  void get_elem( index_t k , std::vector<const real_t*>& X ) const;

  bool has( index_t k , index_t idx ) const;
  bool has( index_t p ) const;

  bool closed() const { return closed_; }
  void set_closed( bool x ) { closed_ = x; }
  void close();

  index_t cardinality( const index_t* v , index_t nv ) const;

  void intersect( const std::vector<index_t>& facet , std::vector<index_t>& elems ) const;

  void all_with( const std::vector<index_t>& facet , std::vector<index_t>& elems ) const;

  void get_boundary( Topology<type>& boundary ) const;

  void get_elements( Topology<type>& topology ) const;

  void facet( const index_t k , const index_t j , std::vector<index_t>& f ) const;

  void orient( real_t* p=NULL );
  void orient( index_t* v , const index_t nv , real_t* p=NULL );
  real_t volume() const;

  void get_volumes( std::vector<real_t>& volumes ) const;

  Neighbours<type>& neighbours() { return neighbours_; }
  const Neighbours<type>& neighbours() const { return neighbours_; }

  InverseTopology<type>& inverse() { return inverse_; }
  const InverseTopology<type>& inverse() const { return inverse_; }

  // local operator functions
  void remove_point( index_t k );
  void apply( Cavity<type>& cavity );

private:
  type master_;

  bool closed_;

  Neighbours<type>      neighbours_;
  InverseTopology<type> inverse_;
};

} // luna

#endif
