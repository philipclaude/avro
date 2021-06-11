//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_ADAPT_H_
#define avro_LIB_ADAPTATION_ADAPT_H_

#include "adaptation/primitive.h"

#include <memory>

namespace avro
{

template<typename type> class Topology;
class Mesh;

class AdaptationParameters;
class MetricAttachment;
template<typename type,typename T> class FieldInterpolation;

struct AdaptationProblem
{
  AdaptationProblem( Mesh& _mesh_in ,
    std::vector<symd<real_t>>& _fld,
    AdaptationParameters& _params ,
    Mesh& _mesh_out ,
    FieldInterpolation<Simplex,Metric>* _interpolation=nullptr ) :
      mesh_in(_mesh_in),
      fld(_fld),
      params(_params),
      mesh_out(_mesh_out),
      interpolation(_interpolation)
  {}
  Mesh& mesh_in; // also modified
  std::vector<symd<real_t>>& fld;
  AdaptationParameters& params;
  Mesh& mesh_out;
  FieldInterpolation<Simplex,Metric>* interpolation;
};
template<typename type> int adapt( AdaptationProblem& problem );

template<typename type>
class AdaptThread
{

public:
  AdaptThread( Topology<type>& topology , MetricField<type>& metric , AdaptationParameters& params );

  void adapt();

  void collapse_edges( bool limit_length=false , bool swapout=false );
  void split_edges( real_t lt , bool limit_length=true , bool swapout=false );
  void swap_edges( real_t qt , index_t npass , bool lcheck=false );
  void smooth_points( index_t nb_iter );

  bool check( const std::string& stage_name ) const { return true; }

private:
  bool swap_edge( index_t p , index_t q , real_t Q0 , real_t lmin0 , real_t lmax0 );

  Topology<type>& topology_;
  MetricField<type>& metric_;

  AdaptationParameters& params_;

  Collapse<type> collapser_;
  Insert<type> inserter_;
  Smooth<type> smoother_;
  EdgeSwap<type> edge_swapper_;
};

} // avro

#endif
