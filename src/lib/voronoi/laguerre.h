#ifndef AVRO_LIB_VORONOI_LAGUERRE_H_
#define AVRO_LIB_VORONOI_LAGUERRE_H_

#include "common/types.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/delaunay.h"

namespace avro
{

namespace delaunay
{

class VoronoiDiagram;

class LaguerreDiagram : public Topology<Polytope>
{
public:
  LaguerreDiagram( const Points& sites , const Topology<Simplex>& domain ,
      const std::vector<real_t>& weights = std::vector<real_t>() );

  void set_target_mass( const std::vector<real_t>& mass );
  void initialize();
  void compute();
  void optimize_cvt();
  void optimize_otm();

  real_t eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW , std::vector<real_t>& volumes );
  real_t eval_objective();

  void set_exact( bool x ) { exact_ = x; }

  void set_delaunay( const real_t* x , coord_t dim );
  void set_weights( const real_t* w );
  void set_volumes( std::vector<real_t>& v ) { volume_ = v; }

  void clear_decomposition();

private:
  Points points_;
  Points domain_points_;
  Topology<Simplex> domain_;
  std::vector<real_t> weight_;
  std::vector<real_t> mass_;
  std::vector<index_t> sites_;
  std::vector<real_t> volume_;

  Delaunay delaunay_;

  bool exact_;

  real_t time_voronoi_;
  real_t time_decompose_;
  real_t time_neighbours_;

  std::shared_ptr<VoronoiDiagram> diagram_;

  Topology<Simplex> decomposition_;
  Points decomposition_points_;
};

} // delaunay

} // avro

#endif
