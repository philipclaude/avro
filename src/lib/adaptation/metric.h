#ifndef avro_LIB_ADAPTATION_METRIC_H_
#define avro_LIB_ADAPTATION_METRIC_H_

#include "common/array.h"
#include "common/table.h"
#include "common/types.h"

#include "mesh/field.h"
#include "mesh/points.h"
#include "mesh/search.h"

#include "numerics/linear_algebra.h"
#include "numerics/matrix.h"

namespace avro
{

template<typename type>
void
interp( const std::vector<type>& alpha ,
             const std::vector<numerics::SymMatrixD<type>>& tensors ,
						 numerics::SymMatrixD<type>& T )
{
	avro_assert( alpha.size()==tensors.size() );

	T = 0;
  for (index_t k=0;k<tensors.size();k++)
	{
		T = T + numerics::logm(tensors[k])*alpha[k];
	}
	T = numerics::expm(T);
}

inline numerics::SymMatrixD<real_t>
intersect( const numerics::SymMatrixD<real_t>& x , const numerics::SymMatrixD<real_t>& y )
{
  coord_t n = x.n();

  numerics::MatrixD<real_t> Z( n , n );
  numerics::MatrixD<real_t> Xinv( n , n );
  numerics::VectorD<real_t> lambda( n );

  // compute inverse of current tensor
  numerics::MatrixD<real_t> X(x);
  numerics::MatrixD<real_t> Y(y);
  Xinv = numerics::inverse(X);

  // compute N = T^{-1} * tensors[k]
  Z = Xinv*Y;

  // compute eigenvectors of N
  numerics::SymMatrixD<real_t> T(Z);

  numerics::MatrixD<real_t> E(n,n);
  numerics::VectorD<real_t> P0(n);

  numpack::DLA::EigenSystem(T,P0,E);

  numerics::MatrixD<real_t> P(n,n);
  P = numpack::DLA::diag(P0);

  numerics::MatrixD<real_t> Pinv(n,n);
  Pinv = numerics::inverse(P);

  // measure the length in direction of each eigenvector
  for (coord_t d=0;d<n;d++)
  {
    numerics::VectorD<real_t> e(n);
    e = P.col(d);
    real_t h1 = 1./std::sqrt( e*X*numpack::Transpose(e) );
    real_t h2 = 1./std::sqrt( e*Y*numpack::Transpose(e) );
    real_t h = std::min(h1,h2);

    lambda(d) = 1./(h*h);
  }

  avro_implement; // this is not unit tested

  // construct the matrix from the eigendecomposition
  return Pinv*numpack::DLA::diag(lambda)*numpack::Transpose(Pinv);
}

class Metric : public numerics::SymMatrixD<real_t>
{
public:
  Metric() :
    numerics::SymMatrixD<real_t>(0),
    number_(0),
    log_(0),
    elem_(0),
    sqdet_(-1)
  {}

  Metric( coord_t number ) :
    numerics::SymMatrixD<real_t>(number),
    number_(number),
    log_(number),
    elem_(0),
    sqdet_(-1)
  {}

  void allocate( coord_t n )
  {
    numerics::SymMatrixD<real_t>::allocate(n);
    log_.allocate(n);
    number_ = n;
  }

  numerics::SymMatrixD<real_t>& log(){ return log_; }
  const numerics::SymMatrixD<real_t>& log() const { return log_; }

  void set( const numerics::SymMatrixD<real_t>& m0 )
  {
    avro_assert( number_ == m0.m() && number_ == m0.n() );
    for (coord_t i=0;i<number_;i++)
    for (coord_t j=i;j<number_;j++)
      (*this)(i,j) = m0(i,j);
  }

  real_t sqdet() const { return sqdet_; }
  index_t elem() const { return elem_; }
  void set_elem( index_t elem ) { elem_  = elem; }

  void calculate()
  {
    log_   = numerics::logm(*this);
    sqdet_ = std::sqrt( numerics::determinant(*this) );
  }

private:
  coord_t number_;
  numerics::SymMatrixD<real_t> log_; // logarithm of this metric
  index_t elem_; // element in some mesh containing this metric
  real_t sqdet_; // sqrt of determinant
};

//
// Metric attachment referenced by MetricField for holding metrics at
// dynamic mesh points
//
class MetricAttachment : public Array<Metric>
{
public:
	template<typename Function>
	MetricAttachment( Function& fn , Points& points ) :
    number_(fn.dim()), points_(points)
	{
		for (index_t k=0;k<points_.nb();k++)
    {
      Metric mk(number_);
      if (points.ghost(k))
      {
        for (index_t i=0;i<number_;i++)
          mk(i,i) = 1.;
        mk.calculate();
        Array<Metric>::add( mk );
      }
      else
			{
        mk.set( fn(points_[k]) );
        mk.calculate();
			  Array<Metric>::add( mk );
        mk.dump();
			}
    }
	}

  MetricAttachment( Points& points , const std::vector<numerics::SymMatrixD<real_t>>& metrics );

	const numerics::SymMatrixD<real_t>& log( const index_t k ) const
    { return Array<Metric>::data_[k].log(); }
	real_t sqdet( index_t k ) const
    { return Array<Metric>::data_[k].sqdet(); }

  template<typename type> void set_cells( const Topology<type>& topology );

  void reset( MetricAttachment& fld );
  void add( numerics::SymMatrixD<real_t>& tensor, index_t elem );

	void assign( index_t p , const numerics::SymMatrixD<real_t>& M , index_t elem );
  void remove( index_t k , bool recheck=true );
  bool check() const;

  template<typename type>
  void limit( const Topology<type>& topology , real_t href );

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
class MetricField : public Field<type,Metric>
{
public:
  MetricField( Topology<type>& topology , MetricAttachment& fld );

	numerics::SymMatrixD<real_t>& operator() ( const Points& x , index_t v );

	real_t length( index_t n0 , index_t n1 ) const;
	real_t length( const Points& v , index_t n0 , index_t n1 ) const
		{ return length(n0,n1); }
  void lengths( const Topology<type>& topology , std::vector<real_t>& lengths ) const;

	real_t volume( const Topology<type>& t , const index_t k );
  real_t volume( const Topology<type>& t );

	real_t quality( const Topology<type>& topology , index_t k );
  int  find( index_t n0 , index_t n1 , real_t*  x );
  void interpolate( real_t*  x , index_t elem , numerics::SymMatrixD<real_t>& tensor , bool STFU=false );
  void add( index_t n0 , index_t n1 , real_t*  x );
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

private:

  // static data
  const Topology<type>& topology_;
  const coord_t number_;

  // the dynamic field that changes when vertices are changed
	MetricAttachment& attachment_;

  ElementSearch<type> searcher_;

  real_t normalization_;
};

template<typename type>
real_t
worst_quality( const Topology<type>& topology , MetricField<type>& metric )
{
  std::vector<real_t> quality(topology.nb());
  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k))
      quality[k] = 1e20;
    else
      quality[k] = metric.quality( topology , k );
  }
  return *std::min_element(quality.begin(),quality.end());
}

} // avro

#endif
