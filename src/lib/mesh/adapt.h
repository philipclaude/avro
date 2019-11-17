#ifndef LUNA_LIB_MESH_ADAPT_H_
#define LUNA_LIB_MESH_ADAPT_H_

namespace luna
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

} // luna

#endif
