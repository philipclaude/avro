#ifndef avro_LIB_LIBRARY_METRIC_H_
#define avro_LIB_LIBRARY_METRIC_H_

#include "common/types.h"

#include "numerics/matrix.h"

namespace avro
{

namespace library
{

class MetricField_Analytic
{
public:
  MetricField_Analytic( coord_t dim ) :
    dim_(dim)
  {}
  virtual ~MetricField_Analytic() {}

  virtual numerics::SymMatrixD<real_t> operator()( const real_t* x ) const = 0;

protected:
  coord_t dim_;
};

class MetricField_Uniform : public MetricField_Analytic
{
public:
  MetricField_Uniform( coord_t dim , real_t h ) :
    MetricField_Analytic(dim),
    h_(dim,h)
  {}

  MetricField_Uniform( coord_t dim , const std::vector<real_t>& h ) :
    MetricField_Analytic(dim),
    h_(h)
  {}

  coord_t dim() const { return dim_; }

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(dim_);

    for (coord_t d=0;d<dim_;d++)
      m(d,d) = 1./(h_[d]*h_[d]);
    return m;
  }

private:
  std::vector<real_t> h_;
};

class MetricField_UGAWG_Linear : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Linear() :
    MetricField_Analytic(3)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    numerics::SymMatrixD<real_t> m(dim_);

    real_t hx = 0.1;
    real_t hy = 0.1;
    real_t h0 = 1e-3;
    real_t hz = h0 +2.*(0.1 -h0)*fabs( x[2] -0.5 );

    m(0,0) = 1./(hx*hx);
    m(0,1) = 0.;
    m(0,2) = 0.;
    m(1,1) = 1./(hy*hy);
    m(1,2) = 0.;
    m(2,2) = 1./(hz*hz);
    return m;
  }

};

} // library

} // avro

#endif
