#ifndef URSA_LIB_COMMON_TREE_H_
#define URSA_LIB_COMMON_TREE_H_

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

  void children( std::vector<Child_ptr>& c ) const;
  void children( std::vector<T*>& c ) const;

  T* child_ptr( index_t k ) { return child_[k].get(); }
  const T* child_ptr( index_t k ) const { return child_[k].get(); }

  T* child_smptr( index_t k ) { return child_[k]; }
  const T* child_smptr( index_t k ) const { return child_[k]; }

  T& child( index_t k ) { return *child_[k].get(); }
  const T& child( index_t k ) const { return *child_[k].get(); }

private:

  std::vector<Child_ptr> child_;
};


} // ursa

#endif
