#ifndef avro_LIB_ADAPTATION_PARALLEL_H_
#define avro_LIB_ADAPTATION_PARALLEL_H_

#include "common/types.h"

#include "mesh/points.h"

#include "numerics/matrix.h"

#include <map>
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

template<typename type>
class AdaptationInterface : public Topology<type>
{
public:
  AdaptationInterface( coord_t dim , coord_t udim , coord_t number ) :
    Topology<type>(points_,number),
    points_(dim,udim)
  {}

  index_t add_point( index_t p , Points& points , index_t offset=0 )
  {
    index_t idx;
    if (global2local_.find(p+offset)==global2local_.end())
    {
      idx = points_.nb();
      global2local_.insert( { p+offset,idx } );
      local2global_.push_back( p+offset );
      points_.create( points[p] );
      points_.set_entity( idx , points.entity(p) );
      points_.set_param( idx , points.u(p) );
    }
    else
      idx = global2local_.at(p+offset);
    return idx;
  }

private:
  Points points_;
  std::vector<index_t> local2global_;
  std::map<index_t,index_t> global2local_;
};

} // avro

#endif
