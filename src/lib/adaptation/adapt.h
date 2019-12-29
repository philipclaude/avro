#ifndef luma_LIB_ADAPTATION_ADAPT_H_
#define luma_LIB_ADAPTATION_ADAPT_H_

#include "adaptation/primitive.h"

#include <memory>

namespace luma
{

template<typename type> class Topology;
class Mesh;

class AdaptationParameters;

struct AdaptationProblem
{
  Mesh& mesh_in; // also modified
  std::vector<numerics::SymMatrixD<real_t>>& fld;
  //VertexField<numerics::SPDT<real>>& fld;
  AdaptationParameters& params;
  Mesh& mesh_out;
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

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( Topology<type>& topology ); // serial adaptation
  AdaptationManager( Mesh& mesh ); // parallel adaptatoin

private:
  std::vector<AdaptThread<type>> thread_;

};

} // luma

#endif
