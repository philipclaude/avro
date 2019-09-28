#include "tree.h"

namespace ursa
{

template<typename Node_t>
void
Tree<Node_t>::addChild( Node_ptr node )
{
  child_.push_back( node );
}

} // ursa
