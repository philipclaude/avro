#ifndef URSA_LIB_COMMON_TREE_H_
#define URSA_LIB_COMMON_TREE_H_

#include "common/error.h"
#include "common/types.h"

#include <memory>
#include <vector>

namespace ursa
{

template<typename Node_t>
class Tree
{
public:
  typedef std::shared_ptr<Node_t> Node_ptr;

  index_t nb_children() const { return child_.size(); }

  void addChild( Node_ptr c );

  void children( std::vector<Node_ptr>& c ) const;
  void children( std::vector<Node_t*>& c ) const;

  Node_t* child_ptr( index_t k ) { return child_[k].get(); }
  const Node_t* child_ptr( index_t k ) const { return child_[k].get(); }

  Node_ptr child_smptr( index_t k )
    { ursa_assert(k<nb_children()); return child_[k]; }
  const Node_ptr child_smptr( index_t k ) const
    { ursa_assert(k<nb_children()); return child_[k]; }

  Node_t& child( index_t k )
    { ursa_assert(k<nb_children()); return *child_[k].get(); }
  const Node_t& child( index_t k ) const
    { ursa_assert(k<nb_children()); return *child_[k].get(); }

private:

  std::vector<Node_ptr> child_;
};


} // ursa

#endif
