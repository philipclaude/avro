//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_SANDBOX_SDOT_MEASURES_H_
#define AVRO_SANDBOX_SDOT_MEASURES_H_

#include "numerics/linear_algebra.h"
#include "numerics/sym.h"
#include "numerics/vec.h"

#include "voronoi/optimal_transport.h"

namespace avro
{
using namespace voronoi;

class DensityMeasure_Sin : public voronoi::DensityMeasure {
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    return 1e1*( 1 + sin(2*M_PI*x[0])*sin(2*M_PI*x[1]) );
  }
};

class DensityMeasure_Quadratic : public voronoi::DensityMeasure {
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    return .1 + 100*x[0]*x[0];
  }
};

class DensityMeasure_Sphere : public voronoi::DensityMeasure {
public:
  DensityMeasure_Sphere( coord_t dim ) :
    dim_(dim)
  {}
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    real_t rho = 1.0;
    for (coord_t d = 0; d < dim_; d++)
      rho += 100.0*(x[d] - 0.5)*(x[d] - 0.5);
    return rho;
  }
private:
  coord_t dim_;
};

class DensityMeasure_Cone : public voronoi::DensityMeasure {
public:
  DensityMeasure_Cone( coord_t dim ) :
    dim_(dim)
  {}
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    real_t r2 = 0.0;
    for (coord_t d = 0; d < dim_-1; d++)
      r2 += x[d]*x[d];
    real_t r = std::sqrt(r2);
    real_t t = x[dim_-1];
    real_t a = 1.0;
    real_t b = r0 - r1;
    real_t c = -r0;
    real_t d = std::abs( a*r + b*t + c ) / std::sqrt( a*a + b*b );
    real_t rho = 100.0/(d*d + 1e-3);
    return rho;
  }
private:
  coord_t dim_;
  const real_t r0 = 0.4;
  const real_t r1 = 0.7;
};


class DensityMeasure_Shock : public voronoi::DensityMeasure {
public:
  DensityMeasure_Shock( coord_t dim ) :
    dim_(dim) {
    k0_ = 10.0;
    k1_ = 1000.0;
    vs_ = 0.7;
    r0_ = 0.4;
    alpha_ = 1e-3;
  }

  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    real_t rt = r0_ + vs_*x[dim_-1];
    real_t x2 = 0.0;
    for (coord_t d = 0; d < dim_-1; d++) {
      x2 += x[d]*x[d];
    }
    x2 = std::sqrt(x2);
    return 1e0 + k0_*std::exp( -alpha_*x[dim_-1] )*std::exp( -k1_*(rt-x2)*(rt-x2) );
  }
private:
  coord_t dim_;
  real_t k0_,k1_,vs_,r0_,alpha_;

};

class DensityMeasure_Gaussian : public voronoi::DensityMeasure {
public:
  DensityMeasure_Gaussian( const vecd<real_t>& mu , const symd<real_t>& sigma ) :
    mu_(mu),
    sigma_(sigma),
    sigma_inv_(mu.m()),
    dim_(mu.m()),
    X_(dim_) {
    detS_ = numerics::det(sigma_);
    sigma_inv_ = numerics::inverse(sigma_);
    avro_assert( dim_ == sigma_.m() );
    avro_assert( dim_ == sigma_.n() );
  }

  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const {
    vecd<real_t> X(dim_);
    for (coord_t d = 0; d < dim_; d++)
      X[d] = x[d] - mu_[d];
    real_t A = 1./std::sqrt( std::pow(2*M_PI,dim_)*detS_);
    real_t exp = numerics::quadratic_form(sigma_inv_,X);
    return A*std::exp(-0.5*exp);
  }

private:
  vecd<real_t> mu_;
  symd<real_t> sigma_;
  symd<real_t> sigma_inv_;
  coord_t dim_;
  vecd<real_t> X_;
  real_t detS_;
};

}

#endif
