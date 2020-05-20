#include "tree.h"

namespace avro
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
Tree<Node_t>::add_parent( Node_t* parent )
{
  parents_.push_back(parent);

  // also add to the children!
  for (index_t k=0;k<nb_children();k++)
    child(k).add_parent(parent);
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

template<typename Node_t>
bool
Tree<Node_t>::has_parent( const Node_t* P ) const
{
  // the children should also have this parent!
  // check them first
  for (index_t k=0;k<nb_children();k++)
  {
    if (!child(k).has_parent(P))
      return false;
  }

  // all the children are deemed to have this parent
  // now check this entity
  for (index_t k=0;k<parents_.size();k++)
  {
    if (parents_[k]==P)
      return true;
  }
  return false;
}

template<typename Node_t>
bool
Tree<Node_t>::has_child( const Node_t* c ) const
{
  for (index_t k=0;k<nb_children();k++)
    if (child_ptr(k)==c) return true;
  return false;
}

template<typename Node_t>
void
Tree<Node_t>::get_children( std::vector<Node_t*>& children )
{
  for (index_t k=0;k<nb_children();k++)
  {
    children.push_back(child_ptr(k));
    child(k).get_children(children);
  }
  uniquify(children);
}

template<typename Node_t>
void
Tree<Node_t>::get_children( std::vector<const Node_t*>& children ) const
{
  for (index_t k=0;k<nb_children();k++)
  {
    children.push_back(child_ptr(k));
    child(k).get_children(children);
  }
  uniquify(children);
}

template<typename Node_t>
void
Tree<Node_t>::get_adjacency( numerics::MatrixD<int>& A ) const
{
  std::vector<const Node_t*> children;
  children.push_back(derived());
  get_children(children);
  A.resize( children.size() , children.size() );
  A = 0;
  for (index_t i=0;i<children.size();i++)
  for (index_t j=i+1;j<children.size();j++)
  {
    if (children[i]->has_child(children[j]))
      A(i,j) = 1;
  }
}

template<typename Node_t>
void
Tree<Node_t>::print( index_t level ) const
{
  for (index_t i=0;i<level;i++) printf("    ");
  derived()->print_header();
  for (index_t k=0;k<nb_children();k++)
    child(k).Tree<Node_t>::print(level+1);
}

template<typename Node_t>
void
Tree<Node_t>::print_header() const
{
  printf("default print!!\n");
}

} // avro
