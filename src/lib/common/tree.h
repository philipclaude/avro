#ifndef luma_LIB_COMMON_TREE_H_
#define luma_LIB_COMMON_TREE_H_

#include "common/error.h"
#include "common/types.h"

#include <memory>
#include <vector>

namespace luma
{

class TreeBase
{

};

template<typename Node_t>
class Tree : public TreeBase
{
public:
  typedef std::shared_ptr<Node_t> Node_ptr;

  void copy( const TreeBase& tree );

  index_t nb_children() const { return child_.size(); }

  void add_child( Node_ptr c );

  void children( std::vector<Node_ptr>& c ) const;
  void children( std::vector<Node_t*>& c ) const;

  bool above( const Node_t* node ) const;

  Node_t* child_ptr( index_t k ) { return child_[k].get(); }
  const Node_t* child_ptr( index_t k ) const { return child_[k].get(); }

  Node_ptr child_smptr( index_t k )
    { luma_assert(k<nb_children()); return child_[k]; }
  const Node_ptr child_smptr( index_t k ) const
    { luma_assert(k<nb_children()); return child_[k]; }

  Node_t& child( index_t k )
    { luma_assert(k<nb_children()); return *child_[k].get(); }
  const Node_t& child( index_t k ) const
    { luma_assert(k<nb_children()); return *child_[k].get(); }

  void build_parents();
  void add_parent( Node_t* parent );
  bool has_parent( const Node_t* P ) const;

  void get_children( std::vector<Node_t*>& children );

protected:
  std::vector<Node_ptr> child_;
  std::vector<Node_t*> parents_; // list of all parents owning this

  Node_t* derived() { return static_cast<Node_t*>(this); }
};


} // luma

#endif
