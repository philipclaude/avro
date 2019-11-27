#ifndef LUNA_LIB_ADAPTATION_ADAPT_H_
#define LUNA_LIB_ADAPTATION_ADAPT_H_

namespace luna
{

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

} // luna

#endif
