#ifndef avro_LIB_ADAPTATION_PARALLEL_H_
#define avro_LIB_ADAPTATION_PARALLEL_H_

#include "common/types.h"

#include "numerics/matrix.h"

#include <vector>

namespace avro
{

template<typename type> class Topology;

typedef numerics::SymMatrixD<real_t> VertexMetric;
class AdaptationParameters;

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( Topology<type>& topology , std::vector<VertexMetric>& metrics , AdaptationParameters& params );

  void adapt();

  void repartition();
  void partition();

private:
  Topology<type>& topology_;
  std::vector<VertexMetric> metrics_;
  AdaptationParameters& params_;
};


template<typename type>
int adaptp( Topology<type>& topology , const std::vector<VertexMetric>& metrics , AdaptationParameters& params , Topology<type>& topology_out , index_t level=0 );

} // avro

#endif
