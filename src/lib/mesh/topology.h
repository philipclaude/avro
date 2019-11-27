#ifndef LUNA_MESH_TOPOLOGY_H_
#define LUNA_MESH_TOPOLOGY_H_

#include "common/array.h"
#include "common/json.h"
#include "common/tree.h"
#include "common/types.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/field.h"

#include <string>
#include <vector>

namespace luna
{

class Points;
class ClippingPlane;

class TopologyHolder : public Array<index_t>
{
public:
  TopologyHolder( Points& vertices , const coord_t number , ArrayLayoutCategory category ) :
    Array<index_t>(category,number+1),
    points_(vertices),
    number_(number),
    fields_(*this)
  {}

  virtual ~TopologyHolder() {}

  Points& points() const { return points_; }

  coord_t number() const { return number_; }
  coord_t order() const { return order_; }

  // virtual functions for leaf
  void copy( TopologyHolder& topology );

  virtual void getPoints( std::vector<index_t>& p ) const = 0;
  virtual void getEdges( std::vector<index_t>& e ) const = 0;
  virtual void getTriangles( std::vector<index_t>& t ) const = 0;

  const std::string& name() const { return name_; }
  void setName( const std::string& _name ) { name_ = _name; }

  void setDummy( bool x ) { dummy_ = x; }
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
};

#if 1

template<typename type>
class Topology : public Tree<Topology<type>>, public TopologyHolder
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

  bool ghost( index_t k ) const { return false; }

  void getPoints( std::vector<index_t>& p ) const {}
  void getEdges( std::vector<index_t>& e ) const;
  void getTriangles( std::vector<index_t>& t ) const;

private:
  type master_;
};

#else

template<typename Shape> class Topology;

template<typename type>
class TopologyBase : public Tree<Topology<type>>, public TopologyHolder
{

public:
  typedef Topology<type> Topology_t;
  typedef std::shared_ptr<Topology_t> Topology_ptr;

  using Tree<Topology_t>::nb_children;

  TopologyBase( Points& _vertices , const coord_t _number );
  TopologyBase( Points& _vertices , const coord_t _number , const coord_t _order );
  TopologyBase( Points& _vertices , const json& J );

  Topology_t& topology( index_t k ) { return Tree<Topology_t>::child(k); }

  type& master() { return master_; }
  const type& master() const { return master_; }

  coord_t order() const { return master_.order(); }

  void do_something() {}

  bool ghost( index_t k ) const { return false; }

  void getPoints( std::vector<index_t>& p ) const {}
  void getEdges( std::vector<index_t>& e ) const;
  void getTriangles( std::vector<index_t>& t ) const;

  type master_;
};


template<>
class Topology<Simplex> : public TopologyBase<Simplex>
{
public:
  Topology( Points& vertices , coord_t order );
  Topology( Points& vertices , const Topology<Simplex>& lagrange );
  Topology( Points& vertices , const Topology<Simplex>& lagrange , coord_t order );

  void convert( const Topology<Simplex>& linear );
};

template<>
class Topology<Polytope> : public TopologyBase<Polytope>
{
public:

private:
  Data<int> incidence_;
};

#endif

} // luna

#endif
