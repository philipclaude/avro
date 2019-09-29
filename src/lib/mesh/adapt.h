#ifndef URSA_LIB_MESH_ADAPT_H_
#define URSA_LIB_MESH_ADAPT_H_

namespace ursa
{

class

template<typename type>
class AdaptThread
{
public:
  AdaptThread( Topology<type>& topology );

  void insertions();
  void collapses();
  void edgeswaps();

  Topology<type>& topology_;
};

} // ursa

#endif
