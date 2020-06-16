//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_IMPLIED_METRIC_H_
#define avro_LIB_ADAPTATION_IMPLIED_METRIC_H_

#include "common/array.h"
#include "common/types.h"

#include "mesh/field.h"

#include "numerics/matrix.h"

namespace avro
{

class Points;
template <typename type> class Topology;

template<typename type>
class ElementImpliedMetric : public numerics::SymMatrixD<real_t>
{
public:
  ElementImpliedMetric( const type&  );

  void compute( const std::vector<const real_t*>& xk );
  void compute( const Points& points , const index_t* v , index_t nv );
  void inverse( const Points& points , const index_t* v , index_t nv  );
  real_t determinant( const std::vector<const real_t*>& xk );
  real_t determinant( const Points& points , const index_t* v , const index_t nv );

private:
  const type& element_;
  numerics::MatrixD<real_t> J_;   // jacobian from physical elem to equilateral
  numerics::MatrixD<real_t> J0_;  // jacobian from physical elem to right-angled
  numerics::MatrixD<real_t> Jeq_; // jacobian from right-angled elem to equilateral
  numerics::MatrixD<real_t> M_;
  real_t detJeq_;
};

template<typename type>
class MeshImpliedMetric : public Array<numerics::SymMatrixD<real_t>>
{
public:
  MeshImpliedMetric( const Topology<type>& topology );

  void initialize();
  void optimize();

  template<int DIM> real_t cost( const std::vector<numerics::SymMatrixD<real_t>>& m , std::vector<numerics::SymMatrixD<real_t>>& dc_dS , real_t& complexity0 ) const;
  template<int DIM> real_t deviation( const std::vector<numerics::SymMatrixD<real_t>>& S , std::vector<numerics::SymMatrixD<real_t>>& df_dS ) const;

  const Topology<type>& topology() const;

private:
  const Topology<type>& topology_;
  std::vector<numerics::SymMatrixD<real_t>> nodalMetricSqrt_;
  std::vector<real_t> nodalMetricSqrtDet_;
  std::vector<index_t> edges_;
  std::vector< std::vector<index_t> > v2e_;
};

} // avro

#endif
