#ifndef LUNA_LIB_ADAPTATION_ADAPT_H_
#define LUNA_LIB_ADAPTATION_ADAPT_H_

#include <memory>

namespace luna
{

template<typename type> class Topology;
template<typename type> class Mesh;

template<typename type>
class AdaptThread
{

public:
  AdaptThread( Topology<type>& topology );

  void adapt();

  void collapse_edges();
  void split_edges();
  void swap_edges();
  void smooth_points();

private:
  Topology<type>& topology_;
};

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( Topology<type>& topology ); // serial adaptation
  AdaptationManager( Mesh<type>& mesh ); // parallel adaptatoin

private:
  std::vector<AdaptThread> thread_;

};

} // luna

#endif
