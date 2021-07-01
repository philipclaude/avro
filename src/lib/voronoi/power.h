#ifndef AVRO_LIB_VORONOI_POWER_H_
#define AVRO_LIB_VORONOI_POWER_H_

#include "avro_types.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

class Lite_Sparse_Matrix;

namespace avro
{

namespace delaunay
{

class VoronoiDiagram;
class PowerDiagram;

class PowerFacets : public Topology<Polytope>
{
public:
  PowerFacets( PowerDiagram& diagram );

  void extract();
  void compute_quantities();


  void construct_laplacian();
  void construct_gradient();
  void construct_divergence();

  void print() const;

  const DOF<real_t>& bij() const { return bij_; }
  const std::vector<real_t>& aij() const { return aij_; }
  const std::map<index_t,Bisector>& indices() const { return indices_; }

private:
  PowerDiagram& diagram_;

  std::vector<index_t> cellL_;
  std::vector<int> cellR_;
  std::vector<int> label_;

  // information needed to construct the hessian
  std::vector<real_t> aij_;
  std::vector<real_t> lij_;
  DOF<real_t> bij_;
  std::map<index_t,Bisector> indices_;

  std::set<int> bisectors_;
  std::map<Bisector,index_t> bisector2facet_;
  std::map<SymbolicVertex,index_t> symbolic_;
  std::map<index_t,index_t> unique_vertex_;
  std::vector<std::set<int>> neighbours_;
};

class PowerDiagram : public Topology<Polytope>
{
public:
  PowerDiagram( const Points& sites , const Topology<Simplex>& domain ,
      const std::vector<real_t>& weights = std::vector<real_t>() );

  void set_target_mass( const std::vector<real_t>& mass );
  void compute();

  void optimize_cvt();
  void optimize_otm();

  real_t eval_objective( std::vector<real_t>& dE_dZ , std::vector<real_t>& dE_dW , std::vector<real_t>& volumes );

  void set_exact( bool x ) { exact_ = x; }

  void set_delaunay( const real_t* x , coord_t dim );
  void set_weights( const real_t* w );
  void set_volumes( std::vector<real_t>& v ) { volume_ = v; }

  void clear_decomposition();

  const VoronoiDiagram& diagram() const { return *diagram_.get(); }

  index_t site( index_t k ) const { avro_assert( k < nb() && k < sites_.size() ); return sites_[k]; }
  const Delaunay& delaunay() const { return delaunay_; }

  void get_hessian( Lite_Sparse_Matrix& h ) const;

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

  PowerFacets facets_;


  std::vector<double> diag_;
  std::vector<double> values_;
  std::vector<int> rowind_;
  std::vector<int> colptr_;
};

} // delaunay

} // avro

#endif
