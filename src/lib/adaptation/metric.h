#ifndef LUNA_LIB_ADAPTATION_METRIC_H_
#define LUNA_LIB_ADAPTATION_METRIC_H_

#include "common/table.h"
#include "common/types.h"

#include "mesh/field.h"

namespace luna
{

#if 0

class Metric : public MatrixSymD<real_t>
{
public:
  MatrixSymD<real_t>& log(){ return log_; }
  const MatrixSymD<real_t>& log() const { return log_; }

  real_t sqdet() const { return sqdet_; }
  index_t elem() const { return elem_; }
  void set_elem( index_t elem ) { elem_  = elem; }

  void calculate()
  {
    log_   = this->log();
    sqdet_ = std::sqrt( numerics::det(*this) );
  }

private:
  index_t elem_;           // element in some mesh containing this metric
  MatrixSymD<real_t> log_; //
  real_t sqdet_;           // sqrt of determinant
};

//
// DiscreteField used by DiscreteMetric for holding metrics
//
class MetricAttachment : public Array<Metric>
{
public:
	template<typename Function>
	MetricAttachment( Function& fn , Points& points ) :
    number_(fn.number()), points_(points)
	{
		for (index_t k=0;k<points_.nb();k++)
    {
      Metric mk(number_);
      if (points.ghost(k))
      {
        for (index_t i=0;i<n_;i++)
          mk(i,i) = 1.;
        mk.calculate();
        Array<Metric>::add( mk );
      }
      else
			{
        mk = fn(points_[k]);
        fn.calculate();
			  Array<Metric>::add( mk );
			}
    }
	}

  //MetricAttachment( VertexField<SPDT<real_t>>& fld , Points& v );
  MetricAttachment( Points& points );

	const MatrixSymD<real_t>& log( const index_t k ) const
    { return metric_[k].log(); }
	real_t sqdet( index_t k ) const
    { return metric_[k].sqdet(); }

  template<typename type> void set_cells( Topology<type>& topology );

  //void reset( DiscreteField& fld );
  //void set( VertexField<SPDT<real_t>>& fld );

  void add( Metric& m , index_t elem );

	void assign( index_t p , const Metric& M , index_t elem );
  void remove( index_t k , bool recheck=true );
  bool check() const;

  Points& points() { return points_; }

  std::vector<index_t>& elems() { return elems_; }
  index_t& elem( const index_t k ) { return elem_[k]; }

  index_t nb_elem() const { return elem_.size(); }

  void to_json( json& J ) const;
	void from_json( json& J ) const;

	void to_solb( const std::string& filename ) const;
	void from_solb( const std::string& filename );

private:
  coord_t n_;
  Vertices& vertices_;
};

//
// discrete metric field
//
template<typename type>
class MetricField : public Field<type,Metric>>
{
public:
  MetricField( coord_t n , Topology<type>& topology , MetricAttachment& fld );
	MetricField( Topology<type>& topology , MetricAttachment& fld );

	Metric<real_t>& operator() ( const Vertices& x , index_t v );
	real length( index_t n0 , index_t n1 );

	real length( const Vertices& v , index_t n0 , index_t n1 )
		{ return length(n0,n1); }

	real volume( const Topology<type>& t , const index_t k );
  real volume( const Topology<type>& t );
  real volume0();

	real quality( const Topology<type>& topology , index_t k );
  void initializeCells();
  int  find( index_t n0 , index_t n1 , real* x );
  void interpolate( real* x , index_t elem , SPDT<real_t>& tensor , bool STFU=false );
  void add( index_t n0 , index_t n1 , real* x );
  void recompute( index_t p , real* x , const std::vector<index_t>& N );
  bool recompute( index_t p , real* x );
	void assign( index_t p , const SPDT<real_t>& M0 , index_t elem0 )
		{ field_.assign(p,M0,elem0); }
	index_t elementContaining( index_t p )
		{ return field_.cell(p); }
  void remove( index_t k );

  template<typename Function>
  void reevaluate( Topology<type>& topology , AnalyticMetric<Function>& analytic );

  void reset( MetricAttachment& fld );

  MetricAttachment& field() { return field_; }
  Topology<type>& topology() { return background_.topology(); }
  ElementSearch<type>& searcher() { return searcher_; }

  bool checkCells();

  // ensures shape = type
  template<typename shape>
  bool check( Topology<shape>& topology );

private:

  // static data
  const Topology<type>& topology_;

  // the dynamic field that changes when vertices are changed
	MetricAttachment& field_;

  ElementSearch<type> searcher_;
};

#endif

} // luna

#endif
