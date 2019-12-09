#include "tree.h"

namespace luna
{

template<typename Node_t>
void
Tree<Node_t>::add_child( Node_ptr node )
{
  child_.push_back( node );
}

template<typename Node_t>
bool
Tree<Node_t>::above( const Node_t* node ) const
{
  for (index_t k=0;k<nb_children();k++)
  {
    const Node_t* ck = child_ptr(k);
    if ( node==ck ) return true;
    if ( ck->above(node) ) return true;
  }
  return false;
}

template<typename Node_t>
void
Tree<Node_t>::build_parents()
{
  // go through the children
  for (index_t k=0;k<nb_children();k++)
  {
    child(k).add_parent(derived());
    child(k).build_parents();
  }
  uniquify(parents_);
}

} // luna
