//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_METRIC_H_
#define avro_LIB_ADAPTATION_METRIC_H_

#include "common/array.h"
#include "common/table.h"
#include "avro_types.h"

#include "mesh/field.h"
#include "mesh/interpolation.h"
#include "mesh/points.h"
#include "mesh/search.h"

#include "numerics/linear_algebra.h"
#include "numerics/sym.h"
#include "numerics/vec.h"

namespace avro
{

template<typename type> class MetricEvaluator;

class Metric : public symd<real_t> {
public:
  Metric() :
    symd<real_t>(0),
    number_(0),
    log_(0),
    elem_(0),
    sqdet_(-1)
  {}

  Metric( coord_t number ) :
    symd<real_t>(number),
    number_(number),
    log_(number),
    elem_(0),
    sqdet_(-1)
  {}

  Metric( const symd<real_t>& A ) :
    Metric(A.m()) {
    set(A);
    calculate();
  }

  void allocate( coord_t n ) {
    symd<real_t>::allocate(n);
    log_.allocate(n);
    number_ = n;
  }

  symd<real_t>& log() { return log_; }
  const symd<real_t>& log() const { return log_; }

  void set( const symd<real_t>& m0 ) {
    avro_assert( number_ == m0.m() && number_ == m0.n() );
    for (coord_t i = 0; i < number_; i++)
    for (coord_t j = i; j < number_; j++)
      (*this)(i,j) = m0(i,j);
  }

  real_t sqdet() const { return sqdet_; }
  index_t elem() const { return elem_; }
  void set_elem( index_t elem ) { elem_  = elem; }

  void calculate() {
    log_   = numerics::logm(*this);
    real_t d = numerics::det(*this);
    avro_assert_msg( d > 0. , "d = %g\n", d );
    sqdet_ = std::sqrt(d);
  }

private:
  coord_t number_;   // topological number
  symd<real_t> log_; // matrix logarithm of this metric
  index_t elem_;     // element in some mesh containing this metric
  real_t sqdet_;     // sqrt of determinant
};

//
// Metric attachment referenced by MetricField for holding metrics at
// dynamic mesh points
//
class MetricAttachment : public Array<Metric> {
public:
	template<typename Function>
	MetricAttachment( Function& fn , Points& points ) :
    number_(fn.dim()), points_(points) {
		for (index_t k=0;k<points_.nb();k++) {
      Metric mk(number_);
      if (points.ghost(k)) {
        for (index_t i=0;i<number_;i++)
          mk(i,i) = 1.;
        mk.calculate();
        Array<Metric>::add( mk );
      }
      else {
        mk.set( fn(points_[k]) );
        mk.calculate();
			  Array<Metric>::add( mk );
			}
    }
	}

  MetricAttachment( Points& points , const std::vector<symd<real_t>>& metrics );
  MetricAttachment( Points& points );

	const symd<real_t>& log( const index_t k ) const
    { return Array<Metric>::data_[k].log(); }
	real_t sqdet( index_t k ) const
    { return Array<Metric>::data_[k].sqdet(); }

  template<typename type> void set_cells( const Topology<type>& topology );

  void reset( MetricAttachment& fld );
  void add( symd<real_t>& tensor, index_t elem );

	void assign( index_t p , const symd<real_t>& M , index_t elem );
  void remove( index_t k , bool recheck=true );
  bool check() const;

  template<typename type>
  void limit( const Topology<type>& topology , real_t href , bool quiet = false );

  Points& points() { return points_; }

  void to_json( json& J ) const;
	void from_json( json& J ) const;

	void to_solb( const std::string& filename ) const;
	void from_solb( const std::string& filename );

private:
  const coord_t number_;
  Points& points_;
};

//
// discrete metric field
//
template<typename type>
class MetricField : public Field<type,Metric> {
public:
  MetricField( Topology<type>& topology , MetricAttachment& fld );
	symd<real_t>& operator() ( const Points& x , index_t v );

	real_t length( index_t n0 , index_t n1 ) const;
	real_t length( const Points& v , index_t n0 , index_t n1 ) const
		{ return length(n0,n1); }
  void lengths( const Topology<type>& topology , std::vector<real_t>& lengths ) const;

	real_t volume( const Topology<type>& t , const index_t k );
  real_t volume( const Topology<type>& t );

	real_t quality( const Topology<type>& topology , index_t k );
  int  find( index_t n0 , index_t n1 , real_t*  x );
  bool add( index_t n0 , index_t n1 , index_t ns , real_t*  x , int idx=-1 );
  bool recompute( index_t p , real_t*  x );
	index_t element_containing( index_t p )
		{ return attachment_[p].elem(); }
  void remove( index_t k );

  void reset( MetricAttachment& fld );

  MetricAttachment& attachment() { return attachment_; }
  const Topology<type>& topology() const { return topology_; }
  ElementSearch<type>& searcher() { return searcher_; }

  bool check_cells();
  bool check( Topology<type>& topology );

  void set_interpolation( FieldInterpolation<type,Metric>* interpolation )
    { interpolation_ = interpolation; }

private:

  const Topology<type>& topology_;
  const coord_t number_;

  // the dynamic field that changes when vertices are changed
	MetricAttachment& attachment_;
  ElementSearch<type> searcher_;
  real_t normalization_;
  FieldInterpolation<type,Metric>* interpolation_;

  mutable std::vector<real_t> edge_;

};

template<typename type>
real_t
worst_quality( const Topology<type>& topology , MetricField<type>& metric , real_t* pquality=nullptr ) {
  std::vector<real_t> quality;//(topology.nb());
  if (pquality == nullptr) {
    quality.resize( topology.nb() );
    pquality = quality.data();
  }
  for (index_t k = 0; k < topology.nb(); k++) {
    if (topology.ghost(k))
      (pquality)[k] = 1e20;
    else
      (pquality)[k] = metric.quality( topology , k );
  }
  return *std::min_element(pquality,pquality + topology.nb());
  //return *std::min_element(quality.begin(),quality.end());
}

} // avro

#endif
