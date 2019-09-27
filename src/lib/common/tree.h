#ifndef URSA_LIB_COMMON_TREE_H_
#define URSA_LIB_COMMON_TREE_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace ursa
{

class TreeNodeBase
{

};

template<typename T>
class Tree
{
public:
  typedef std::shared_ptr<T> Child_ptr;

  index_t nb_children() const { return child_.size(); }

  void addChild( Child_ptr c )
  {
    child_.push_back(c);
  }

  void children( std::vector<Child_ptr>& c ) const;
  void children( std::vector<T*>& c ) const;

  T* child_ptr( index_t k ) { return child_[k].get(); }
  const T* child_ptr( index_t k ) const { return child_[k].get(); }

  T* child_smptr( index_t k ) { return child_[k]; }
  const T* child_smptr( index_t k ) const { return child_[k]; }

  T& child( index_t k ) { ursa_assert(k<nb_children()); return *child_[k].get(); }
  const T& child( index_t k ) const { ursa_assert(k<nb_children()); return *child_[k].get(); }

private:

  std::vector<Child_ptr> child_;
};


} // ursa

#endif
