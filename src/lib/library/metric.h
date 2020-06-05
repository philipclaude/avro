#ifndef avro_LIB_LIBRARY_METRIC_H_
#define avro_LIB_LIBRARY_METRIC_H_

#include "common/error.h"
#include "common/types.h"

#include "mesh/interpolation.h"

#include "numerics/matrix.h"

namespace avro
{

class Metric;

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

  coord_t dim() const { return dim_; }

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

class MetricField_UGAWG_Polar1 : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Polar1() :
    MetricField_Analytic(3)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    real_t r = std::sqrt( x[0]*x[0] +x[1]*x[1] );
    real_t t = atan2( x[1] , x[0] );

    real_t hz = 0.1;
    real_t ht = 0.1;
    real_t h0 = 1e-3;
    real_t hr = h0 +2.*(0.1 -h0)*fabs( r -0.5 );

    numerics::MatrixD<real_t> Q(3,3);
    Q(0,0) = cos(t);
    Q(0,1) = -sin(t);
    Q(0,2) = 0.0;

    Q(1,0) = sin(t);
    Q(1,1) = cos(t);
    Q(1,2) = 0.0;

    Q(2,0) = 0.0;
    Q(2,1) = 0.0;
    Q(2,2) = 1.0;

    numerics::VectorD<real_t> lambda(3);
    lambda[0] = 1./(hr*hr);
    lambda[1] = 1./(ht*ht);
    lambda[2] = 1./(hz*hz);

    numerics::SymMatrixD<real_t> m(3);
    m = Q*numpack::DLA::diag(lambda)*numpack::Transpose(Q);

    return m;
  }
};


class MetricField_UGAWG_Polar2 : public MetricField_Analytic
{
public:
  MetricField_UGAWG_Polar2() :
    MetricField_Analytic(3)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    real_t r = std::sqrt( x[0]*x[0] +x[1]*x[1] );
    real_t t = atan2( x[1] , x[0] );

    real_t hz = 0.1;
    real_t ht = 0.1;
    real_t h0 = 1e-3;
    real_t hr = h0 +2.*(0.1 -h0)*fabs( r -0.5 );
  	real_t d = 10*(0.6 -r);
  	if (d<0.) ht = 0.1;
  	else ht = d/40. +0.1*( 1. -d);

    numerics::MatrixD<real_t> Q(3,3);
    Q(0,0) = cos(t);
    Q(0,1) = -sin(t);
    Q(0,2) = 0.0;

    Q(1,0) = sin(t);
    Q(1,1) = cos(t);
    Q(1,2) = 0.0;

    Q(2,0) = 0.0;
    Q(2,1) = 0.0;
    Q(2,2) = 1.0;

    numerics::VectorD<real_t> lambda(3);
    lambda[0] = 1./(hr*hr);
    lambda[1] = 1./(ht*ht);
    lambda[2] = 1./(hz*hz);

    numerics::SymMatrixD<real_t> m(3);
    m = Q*numpack::DLA::diag(lambda)*numpack::Transpose(Q);

    return m;
  }
};

class MetricField_Tesseract_Linear : public MetricField_Analytic
{
public:
  static constexpr real_t hmin_default = 0.0025;

  MetricField_Tesseract_Linear( real_t hmin = hmin_default ) :
    MetricField_Analytic(3),
    hmin_(hmin)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    real_t hx = 100*hmin_;
    real_t hy = 100*hmin_;
    real_t hz = 100*hmin_;
    real_t h0 = hmin_;
    real_t ht = h0 +2.*(hx -h0)*fabs( x[3] -0.5 );

    numerics::SymMatrixD<real_t> m(4);
    m(0,0) = 1./(hx*hx);
    m(0,1) = 0.;
    m(0,2) = 0.;
    m(0,3) = 0.;
    m(1,1) = 1./(hy*hy);
    m(1,2) = 0.;
    m(1,3) = 0.;
    m(2,2) = 1./(hz*hz);
    m(2,3) = 0.;
    m(3,3) = 1./(ht*ht);
    return m;
  }
private:
  real_t hmin_;
};

class MetricField_Tesseract_Wave : public MetricField_Analytic
{
public:
  MetricField_Tesseract_Wave() :
    MetricField_Analytic(3)
  {}

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const
  {
    real_t eps = 0.001; // offset for singularity at origin
  	real_t X = x[0] +eps;
  	real_t Y = x[1] +eps;
  	real_t Z = x[2] +eps;
  	real_t T = x[3];

    real_t RHO0 = std::sqrt(3*eps*eps); // offset radius from origin

  	// initial and final blast radii
  	real_t r0 = 0.4 +RHO0;
  	real_t rf = 0.8 +RHO0;

  	// hypercone angle
  	real_t alpha = atan2(1.0,(rf-r0));

  	// spherical coordinates
  	real_t RHO = sqrt(X*X+Y*Y+Z*Z);
  	real_t THETA = acos(Z/RHO);
  	real_t PHI = atan2(Y,X);

  	// eigenvectors are normal and tangents to the hypercone
  	numerics::MatrixD<real_t> Q(4,4);
  	Q(0,0) =  sin(alpha)*cos(PHI)*sin(THETA);
  	Q(0,1) =  cos(PHI)*cos(THETA);
  	Q(0,2) = -sin(PHI);
  	Q(0,3) =  cos(alpha)*cos(PHI)*sin(THETA);

  	Q(1,0) =  sin(alpha)*sin(PHI)*sin(THETA);
  	Q(1,1) =  cos(THETA)*sin(PHI);
  	Q(1,2) =  cos(PHI);
  	Q(1,3) =  cos(alpha)*sin(PHI)*sin(THETA);

  	Q(2,0) =  sin(alpha)*cos(THETA);
  	Q(2,1) = -sin(THETA);
  	Q(2,2) = 0.0;
  	Q(2,3) =  cos(alpha)*cos(THETA);

  	Q(3,0) = -cos(alpha);
  	Q(3,1) = 0.0;
  	Q(3,2) = 0.0;
  	Q(3,3) =  sin(alpha);

  	real_t rho0 = r0 +(rf-r0)*T; // blast speed is derivative wrt T
  	real_t h0 = 0.0025;
  	real_t hu = 0.125;
  	real_t hmin = 0.05;
  	real_t delta = 0.1;
  	real_t ht = 0.5;

  	numerics::VectorD<real_t> L(4);
  	real_t hrho = h0 +2*(hu-h0)*fabs(RHO-rho0);
  	real_t htheta = hu;
  	real_t hphi = hu;

  	if (fabs(RHO-rho0)>delta)
  		htheta = hu;
  	else htheta = (hu -hmin)*fabs(RHO-rho0)/delta +hmin;
  	hphi = htheta;

  	L[0] = 1./(hrho*hrho);
  	L[1] = 1./(htheta*htheta);
  	L[2] = 1./(hphi*hphi);
  	L[3] = 1./(ht*ht);

  	return Q*numpack::DLA::diag(L)*numpack::Transpose(Q);
  }
};

template<typename type>
class MetricField_UniformGeometry : public MetricField_Analytic, public FieldInterpolation<type,Metric>
{
public:
  MetricField_UniformGeometry( coord_t dim , real_t hfactor );

  int eval( const Points& points , index_t p , const std::vector<index_t>& guesses , Metric& mp ) override;

  numerics::SymMatrixD<real_t> operator()( const real_t* x ) const override
  {
    avro_assert_not_reached;
  }

  numerics::SymMatrixD<real_t> operator()( const Points& points , index_t p );

private:
  using FieldInterpolation<type,Metric>::analytic_;
  real_t hfactor_;
};

} // library

} // avro

#endif
