#ifndef URSA_MESH_TOPOLOGY_H_
#define URSA_MESH_TOPOLOGY_H_

#include "common/data.h"
#include "common/json.h"
#include "common/tree.h"
#include "common/types.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include <vector>

namespace ursa
{

class Vertices;
class ClippingPlane;

class TopologyHolder : public Data<index_t>
{
public:
  TopologyHolder( Vertices& vertices , const coord_t number ) :
    Data<index_t>(true) , vertices_(vertices) , number_(number)
  {}

  virtual ~TopologyHolder() {}

  Vertices& vertices() const { return vertices_; }

  coord_t number() const { return number_; }
  coord_t order() const { return order_; }

  // virtual functions for leaf
  void copy( TopologyHolder& topology );

  // index/cell retrieval
  std::vector<index_t> get( const index_t k ) const
    { return Data::get(k); }
  index_t operator() ( const index_t k , const index_t j ) const
    { return Data<index_t>::operator()(k,j); }
  index_t& operator() ( const index_t k , const index_t j )
    { return Data<index_t>::operator()(k,j); }

  virtual void getPoints( std::vector<index_t>& p ) const = 0;
  virtual void getEdges( std::vector<index_t>& e ) const = 0;
  virtual void getTriangles( std::vector<index_t>& t ) const = 0;

  const std::string& name() const { return name_; }
  void setName( const std::string& _name ) { name_ = _name; }

  void setDummy( bool x ) { dummy_ = x; }
  bool dummy() const { return dummy_; }

protected:
  Vertices& vertices_;
  coord_t number_;
  coord_t order_;
  bool dummy_;
  std::string name_;
};

template<typename Master_t> class Topology;
template<typename Basis> class Topology<Simplex<Basis>>;

template <typename type>
class TopologyBase : public Tree<Topology<type> >, public TopologyHolder
{

public:
  typedef Topology<type> Topology_t;
  typedef std::shared_ptr<Topology_t> Topology_ptr;

  using Tree<Topology_t>::nb_children;

  TopologyBase( Vertices& _vertices , const coord_t _number );
  TopologyBase( Vertices& _vertices , const coord_t _number , const coord_t _order );
  TopologyBase( Vertices& _vertices , const json& J );

  Topology_t& topology( index_t k ) { return Tree<Topology_t>::child(k); }

  type& master() { return master_; }
  const type& master() const { return master_; }

  coord_t number() const { return master_.number(); }
  coord_t order() const { return master_.order(); }

  void do_something() {}

  void getPoints( std::vector<index_t>& p ) const {}
  void getEdges( std::vector<index_t>& e ) const {}
  void getTriangles( std::vector<index_t>& t ) const {}

  type master_;
};


template<>
class Topology< Simplex<Lagrange> > : public TopologyBase<Simplex<Lagrange>>
{
public:
  typedef Simplex<Lagrange> Master_t;
  using TopologyBase<Master_t>::TopologyBase;
  using TopologyBase<Master_t>::master_;

  Topology( Vertices& vertices , const Topology<Master_t>& linear , coord_t order );

  void convert( const Topology<Master_t>& linear );
};

template<>
class Topology< Simplex<Bezier> > : public TopologyBase<Simplex<Bezier>>
{
public:
  Topology( Vertices& vertices , const Topology<Simplex<Lagrange>>& lagrange );

  void convert();

private:
  const Topology<Simplex<Lagrange>>& lagrange_;
};

template<>
class Topology<Polytope> : public TopologyBase<Polytope>
{
private:
  Data<int> incidence_;
};

}

#endif