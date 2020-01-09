#ifndef avro_LIB_COMMON_TREE_H_
#define avro_LIB_COMMON_TREE_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/matrix.h"

#include <memory>
#include <vector>

namespace avro
{

class TreeBase
{
public:
  /*
  virtual ~TreeBase() {}
  virtual void get_children( std::vector<const TreeBase*>& all ) const = 0;
  virtual index_t nb_children() const = 0;

  template<typename Node_t>
  void construct( std::shared_ptr<Node_t> node ) const
  {
    node = std::make_shared<Node_t>(*this);
  }

private:
  virtual index_t indexof( index_t k , const std::vector<const TreeBase*>& children ) const = 0;
  */
};

template<typename Node_t>
class Tree : public TreeBase
{
public:
  typedef std::shared_ptr<Node_t> Node_ptr;
  ~Tree() {}

  index_t nb_children() const { return child_.size(); }

  template<typename Friend_t> void copy( const Friend_t& tree );

  void add_child( Node_ptr c );

  void children( std::vector<Node_ptr>& c ) const;
  void children( std::vector<Node_t*>& c ) const;

  bool above( const Node_t* node ) const;

  Node_t* child_ptr( index_t k ) { return child_[k].get(); }
  const Node_t* child_ptr( index_t k ) const { return child_[k].get(); }

  Node_ptr child_smptr( index_t k )
    { avro_assert(k<nb_children()); return child_[k]; }
  const Node_ptr child_smptr( index_t k ) const
    { avro_assert(k<nb_children()); return child_[k]; }

  Node_t& child( index_t k )
    { avro_assert(k<nb_children()); return *child_[k].get(); }
  const Node_t& child( index_t k ) const
    { avro_assert(k<nb_children()); return *child_[k].get(); }

  bool has_child( const Node_t* c ) const;

  void build_parents();
  void add_parent( Node_t* parent );
  bool has_parent( const Node_t* P ) const;

  void get_children( std::vector<Node_t*>& children );
  void get_children( std::vector<const Node_t*>& children ) const;
  void get_adjacency( const std::vector<const Node_t*>& children , numerics::MatrixD<index_t>& A ) const;

private:
  index_t indexof( index_t k , const std::vector<const TreeBase*>& children ) const;

protected:
  std::vector<Node_ptr> child_;
  std::vector<Node_t*> parents_; // list of all parents owning this

  Node_t* derived() { return static_cast<Node_t*>(this); }
};

template<typename Node_t>
template<typename Friend_t>
void
Tree<Node_t>::copy( const Friend_t& tree )
{
  std::vector<const Friend_t*> children0;
  children0.push_back( &tree );
  tree.get_children(children0);

  numerics::MatrixD<index_t> A( children0.size() , children0.size() );
  A = index_t(0);
  tree.get_adjacency( children0 , A );

  std::vector<std::shared_ptr<Node_t>> children1( children0.size() );
  for (index_t k=0;k<children0.size();k++)
    children0[k]->construct( children1[k] , *derived() );

  // assign the children
  for (index_t i=0;i<A.m();i++)
  for (index_t j=0;j<A.n();j++)
  {
    if (A(i,j)==0) continue;
    if (i==0)
      add_child( children1[ A(i,j) ] );
    else
      children1[i]->add_child( children1[ A(i,j) ] );
  }
}

} // avro

#endif
