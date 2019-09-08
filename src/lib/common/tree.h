#ifndef URSA_LIB_COMMON_TREE_H_
#define URSA_LIB_COMMON_TREE_H_

#include <memory>
#include <vector>

namespace ursa
{

template<typename T>
class Tree
{
public:

private:
  typedef std::shared_ptr<T> Child_ptr;

  std::vector<Child_ptr> child_;
};


} // ursa

#endif
