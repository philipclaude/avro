#ifndef URSA_LIB_COMMON_GRAPH_H_
#define URSA_LIB_COMMON_GRAPH_H_

namespace ursa
{

template<typename Derived_t>
class Graph
{
public:

  void path( const Node_t& node_i , const Node_t& node_j ) const;

  real_t distance( const Node_t& node_i , const Node_t& node_j ) const;

private:
  Derived_t& derived() { return *static_cast<Derived_t*>(this); }
  const Derived& derived() { return *static_cast<const Derived_t*>(this); }

};

} // ursa

#endif
