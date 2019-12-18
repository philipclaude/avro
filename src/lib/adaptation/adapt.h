#ifndef LUNA_LIB_ADAPTATION_ADAPT_H_
#define LUNA_LIB_ADAPTATION_ADAPT_H_

#include "adaptation/primitive.h"

#include <memory>

namespace luna
{

template<typename type> class Topology;
template<typename type> class Mesh;

class AdaptationParameters;

template<typename type>
class AdaptThread
{

public:
  AdaptThread( Topology<type>& topology , MetricField<type>& metric );

  void adapt();

  void collapse_edges();
  void split_edges();
  void swap_edges( real_t qt , index_t npass , bool lcheck );
  void smooth_points();

private:
  Topology<type>& topology_;
  MetricField<type>& metric_;

  AdaptationParameters& params_;

  Collapse<type> collapser_;
  Insert<type> inserter_;
  Smooth<type> smoother_;
  EdgeSwap<type> edge_swapper_;
  FacetSwap<type> facet_swapper_;
  RidgeSwap<type> ridge_swapper_;
};

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( Topology<type>& topology ); // serial adaptation
  AdaptationManager( Mesh<type>& mesh ); // parallel adaptatoin

private:
  std::vector<AdaptThread<type>> thread_;

};

} // luna

#endif
