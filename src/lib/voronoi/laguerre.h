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
  void compute();
  void optimize();

  real_t eval_objective( std::vector<real_t>& gradient ) const;

  void set_exact( bool x ) { exact_ = x; }

private:
  Points points_;
  Points domain_points_;
  Topology<Simplex> domain_;
  std::vector<real_t> weight_;
  std::vector<real_t> mass_;

  Delaunay delaunay_;

  bool exact_;

};

} // delaunay

} // avro

#endif
