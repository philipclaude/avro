#ifndef avro_LIB_COMMON_TREE_H_
#define avro_LIB_COMMON_TREE_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/matrix.h"

#include <memory>
#include <type_traits>
#include <vector>

namespace avro
{

class TreeBase
{};

template<typename Node_t>
class Tree : public TreeBase
{
public:
  typedef std::shared_ptr<Node_t> Node_ptr;
  //~Tree() {}

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
  void get_adjacency( numerics::MatrixD<int>& A ) const;
  template <typename T> void get_children_typed( std::vector<const T*>& children ) const
  {
    static_assert( std::is_base_of<T,Node_t>::value , "bad types" );
    std::vector<const Node_t*> children0;
    get_children(children0);
    for (index_t k=0;k<children0.size();k++)
      children.push_back( static_cast<const T*>(children0[k]) );
  }

  void print( index_t level=0 ) const;
  void print_header() const;

protected:
  std::vector<Node_ptr> child_;
  std::vector<Node_t*> parents_; // list of all parents owning this

  Node_t* derived() { return static_cast<Node_t*>(this); }
  const Node_t* derived() const { return static_cast<const Node_t*>(this); }
};

template<typename Node_t>
template<typename Friend_t>
void
Tree<Node_t>::copy( const Friend_t& tree )
{
  std::vector<const Friend_t*> children0;
  tree.get_children(children0);
  children0.insert( children0.begin() , &tree );

  numerics::MatrixD<int> A;
  tree.get_adjacency(A);

  A.dump();

  std::vector<std::shared_ptr<Node_t>> children1( children0.size() );
  for (index_t k=0;k<children0.size();k++)
    children0[k]->construct( children1[k] , *derived() );

  // assign the children
  for (index_t i=0;i<index_t(A.m());i++)
  for (index_t j=0;j<index_t(A.n());j++)
  {
    if (A(i,j)==0) continue;
    if (i==0)
      add_child( children1[j] ); // add to root
    else
      children1[i]->add_child( children1[j] ); // add to leaf
  }
}

} // avro

#endif
