#ifndef LUNA_LIB_LIBRARY_METRIC_H_
#define LUNA_LIB_LIBRARY_METRIC_H_

#include "common/types.h"

#include "numerics/matrix.h"

namespace luna
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

} // library

} // luna

#endif
