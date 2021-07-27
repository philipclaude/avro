//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_METRIC_H_
#define avro_LIB_LIBRARY_METRIC_H_

#include "common/error.h"
#include "avro_types.h"

#include "mesh/interpolation.h"

#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/sym.h"

namespace avro
{

class Metric;
class Entity;

namespace library
{

class MetricField_Analytic
{
public:
  MetricField_Analytic( coord_t dim ) :
    dim_(dim),
    scale_(1.0)
  {}
  virtual ~MetricField_Analytic() {}

  virtual symd<real_t> operator()( const real_t* x ) const = 0;

  coord_t dim() const { return dim_; }
  void set_scale( real_t scale ) { scale_ = scale; }
  void scale( matd<real_t>& m ) {
    for (coord_t i = 0; i < dim_; i++)
    for (coord_t j = 0; j < dim_; j++)
      m(i,j) *= scale_;
  }

protected:
  coord_t dim_;
  real_t scale_;
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

  symd<real_t> operator()( const real_t* x ) const
  {
    symd<real_t> m(dim_);

    for (coord_t d=0;d<dim_;d++)
      m(d,d) = 1./(h_[d]*h_[d]);
    return m;
  }

private:
  std::vector<real_t> h_;
};


class MetricField_UGAWG_Linear2d : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Linear2d() :
    MetricField_Analytic(2)
  {}

  symd<real_t> operator()( const real_t* x ) const
  {
    symd<real_t> m(dim_);

    real_t hu = 0.1;
    real_t h0 = hu/100;
    real_t hy = h0 +2.*(hu -h0)*fabs( x[1] -0.5 );
    real_t hx = hu;;

    m(0,0) = 1./(hx*hx);
    m(0,1) = 0.;
    m(1,1) = 1./(hy*hy);
    return m;
  }
};

class MetricField_UGAWG_Linear : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Linear() :
    MetricField_Analytic(3)
  {}

  symd<real_t> operator()( const real_t* x ) const
  {
    symd<real_t> m(dim_);

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

class MetricField_UGAWG_Polar1 : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Polar1() :
    MetricField_Analytic(3)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};


class MetricField_UGAWG_Polar2 : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Polar2() :
    MetricField_Analytic(3)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};

class MetricField_HyperCylinder : public MetricField_Analytic
{
public:
  MetricField_HyperCylinder() :
    MetricField_Analytic(4)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};

class MetricField_Tesseract_Linear : public MetricField_Analytic
{
public:
  static constexpr real_t hmin_default = 0.0025;

  MetricField_Tesseract_Linear( real_t hmin = hmin_default ) :
    MetricField_Analytic(4),
    hmin_(hmin)
  {}

  symd<real_t> operator()( const real_t* x ) const;
private:
  real_t hmin_;
};

class MetricField_Tesseract_RotatingBoundaryLayer : public MetricField_Analytic
{
public:
  static constexpr real_t hmin_default = 0.0025;

  MetricField_Tesseract_RotatingBoundaryLayer() :
    MetricField_Analytic(4)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};

class MetricField_Cube_RotatingBoundaryLayer : public MetricField_Analytic
{
public:
  static constexpr real_t hmin_default = 0.0025;

  MetricField_Cube_RotatingBoundaryLayer() :
    MetricField_Analytic(3)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};


class MetricField_Cube_Wave : public MetricField_Analytic
{
public:
  MetricField_Cube_Wave() :
    MetricField_Analytic(3)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};


class MetricField_Tesseract_Wave : public MetricField_Analytic
{
public:
  MetricField_Tesseract_Wave() :
    MetricField_Analytic(4)
  {}

  symd<real_t> operator()( const real_t* x ) const;
};

template<typename type>
class MetricField_UniformGeometry : public MetricField_Analytic, public FieldInterpolation<type,Metric>
{
public:
  MetricField_UniformGeometry( coord_t dim , real_t hfactor );

  int eval( const Points& points , index_t p , const std::vector<index_t>& guesses , Metric& mp ) override;
  int eval_face( const Points& points , index_t p , Entity* entity , Metric& mp );

  symd<real_t> operator()( const real_t* x ) const override
  {
    avro_assert_not_reached;
  }

  symd<real_t> operator()( const Points& points , index_t p );

private:
  using FieldInterpolation<type,Metric>::analytic_;
  real_t hfactor_;
};

} // library

} // avro

#endif
