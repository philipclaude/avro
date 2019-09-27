#ifndef URSA_MESH_TOPOLOGY_H_
#define URSA_MESH_TOPOLOGY_H_

#include "common/data.h"
#include "common/json.h"
#include "common/tree.h"
#include "common/types.h"

#include "master/master.h"

#include <vector>

namespace ursa
{

class Vertices;
class ClippingPlane;

class TopologyBase : public Data<index_t>, public TreeNodeBase
{
public:
  TopologyBase( Vertices& vertices , const coord_t number ) :
    Data<index_t>(true) , vertices_(vertices) , number_(number)
  {}

  virtual ~TopologyBase() {}

  Vertices& vertices() const { return vertices_; }

  coord_t number() const { return number_; }
  void set_number( const coord_t _number ) { number_ = _number; }

  // virtual functions for leaf
  void copy( TopologyBase& topology );

  // index/cell retrieval
  std::vector<index_t> get( const index_t k ) const
    { return Data::get(k); }
  index_t operator() ( const index_t k , const index_t j ) const
    { return Data<index_t>::operator()(k,j); }
  index_t& operator() ( const index_t k , const index_t j )
    { return Data<index_t>::operator()(k,j); }

  /*virtual void getChildren( std::vector<TopologyBase*>& children ) const = 0;
  virtual void getPoints( std::vector<index_t>& p ) const = 0;
  virtual void getEdges( std::vector<index_t>& e , ClippingPlane* plane ) const = 0;
  virtual void getTriangles( std::vector<index_t>& t , Vertices& v, std::vector<index_t>& parent , ClippingPlane* plane ) const = 0;*/

  const std::string& name() const { return name_; }
  void setName( const std::string& _name ) { name_ = _name; }

  void setDummy( bool x ) { dummy_ = x; }
  bool dummy() const { return dummy_; }

protected:
  Vertices& vertices_;
  coord_t number_;
  bool dummy_;
  std::string name_;
};

template<typename type> class Topology;

template <typename type>
class _Topology : public Tree<Topology<type> >, public TopologyBase
{

public:
  typedef Topology<type> Topology_t;
  typedef std::shared_ptr<Topology_t> Topology_ptr;

  using Tree<Topology_t>::nb_children;

  _Topology( Vertices& _vertices , const coord_t _number );
  _Topology( Vertices& _vertices , const coord_t _number , const coord_t _order );
  _Topology( Vertices& _vertices , const json& J );

  Topology_t& child( index_t k ) { return Tree<Topology_t>::child(k); }

  void do_something() {}

  type master_;
};

template<typename Basis>
class Topology< Simplex<Basis> > : public _Topology< Simplex<Basis> >
{
public:
  Topology( Vertices& vertices );
};

template<>
class Topology<ConvexPolytope> : public _Topology<ConvexPolytope>
{

private:
  Data<int> facets_;
};

}

#endif
