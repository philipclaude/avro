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

  real_t eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW ) const;
  real_t eval_objective() const;

  void set_exact( bool x ) { exact_ = x; }

  void set_delaunay( const real_t* x , coord_t dim );

private:
  Points points_;
  Points domain_points_;
  Topology<Simplex> domain_;
  std::vector<real_t> weight_;
  std::vector<real_t> mass_;
  std::vector<index_t> sites_;

  Delaunay delaunay_;

  bool exact_;

};

} // delaunay

} // avro

#endif
