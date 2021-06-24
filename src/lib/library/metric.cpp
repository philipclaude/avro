//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/metric.h"

#include "geometry/egads/object.h"

#include "library/metric.h"

#include <egads.h>

namespace avro
{

namespace library
{

template<typename type>
MetricField_UniformGeometry<type>::MetricField_UniformGeometry( coord_t dim , real_t h ) :
  MetricField_Analytic(dim),
  FieldInterpolation<type,Metric>(nullptr),
  hfactor_(h)
{
  analytic_ = true;
}

template<typename type>
int
MetricField_UniformGeometry<type>::eval( const Points& points , index_t p , const std::vector<index_t>& guesses , Metric& mp )
{
  symd<real_t> m(dim_);// = uniform_(points[p]);
  Entity* entity = points.entity(p);
  avro_assert_msg( entity!=nullptr , "entity is null for vertex %lu?" , p );

  if (entity->number()==2)
  {
    // evaluate the metric on the face
    eval_face(points,p,entity,mp);
    return 0;
  }

  // interpolate the metrics evaluated on each face
  std::vector<symd<real_t>> metrics;
  for (index_t k=0;k<entity->nb_parents();k++)
  {
    Entity* face = entity->parents(k);
    if (face->number()==2 && face->tessellatable())
    {
      Metric mk(mp.n());
      eval_face( points , p , face , mk );
      metrics.push_back(mk);
    }
  }
  avro_assert( metrics.size()!=0 );
  std::vector<real_t> alpha(metrics.size(),1./metrics.size());
  mp = numerics::interp( alpha , metrics );
  return 0;
}

template<typename type>
int
MetricField_UniformGeometry<type>::eval_face( const Points& points , index_t p , Entity* entity , Metric& mp )
{
  avro_assert( entity->number()==2 );

  symd<real_t> m(dim_);
  real_t area,h;
  real_t range[4];
  int periodic;

  EGADS::Object* eg = (EGADS::Object*) entity;
  #ifndef AVRO_NO_ESP
  EGADS_ENSURE_SUCCESS( EG_getArea( *eg->object() , NULL , &area ) );
  #else
  area = 1.0;
  #endif
  EGADS_ENSURE_SUCCESS( EG_getRange( *eg->object() , range , &periodic ) );

  real_t lu = range[1] - range[0];
  real_t lv = range[3] - range[2];
  real_t area_u = lu*lv;
  h = std::sqrt(area_u);

  h = h*hfactor_;

  m(0,0) = 1./(h*h);
  m(1,1) = 1./(h*h);
  mp.set(m);

  return 0; // returned element should not be used
}

template<typename type>
symd<real_t>
MetricField_UniformGeometry<type>::operator()( const Points& points , index_t p )
{
  Metric m(dim_);
  if (p < points.nb_ghost())
  {
    m = 0.0;
    return m;
  }
  eval( points , p , {} , m );
  return m;
}

symd<real_t>
MetricField_Tesseract_Linear::operator()( const real_t* x ) const {
  real_t hx = 100*hmin_;
  real_t hy = 100*hmin_;
  real_t hz = 100*hmin_;
  real_t h0 = hmin_;
  real_t ht = h0 +2.*(hx -h0)*fabs( x[3] -0.5 );

  symd<real_t> m(4);
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

symd<real_t>
MetricField_Tesseract_Wave::operator()( const real_t* x ) const {
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
  matd<real_t> Q(4,4);
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

  vecd<real_t> L(4);
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

  for (coord_t d = 0; d < 4; d++)
    L[d] *= 2;

  std::pair< vecd<real_t> , matd<real_t> > decomp = {L,Q};
  symd<real_t> m(decomp);
  return m;
}

symd<real_t>
MetricField_UGAWG_Polar1::operator()( const real_t* X ) const {
  real_t eps = 1e-3;
  real_t x[3] = { X[0]+eps, X[1]+eps, X[2] };
  real_t r = std::sqrt( x[0]*x[0] +x[1]*x[1] );
  real_t t = atan2( x[1] , x[0] );

  real_t hz = 0.1;
  real_t ht = 0.1;
  real_t h0 = 1e-3;
  real_t hr = h0 +2.*(0.1 -h0)*fabs( r -0.5 );

  matd<real_t> Q(3,3);
  Q(0,0) = cos(t);
  Q(0,1) = -sin(t);
  Q(0,2) = 0.0;

  Q(1,0) = sin(t);
  Q(1,1) = cos(t);
  Q(1,2) = 0.0;

  Q(2,0) = 0.0;
  Q(2,1) = 0.0;
  Q(2,2) = 1.0;

  vecd<real_t> lambda(3);
  lambda[0] = 1./(hr*hr);
  lambda[1] = 1./(ht*ht);
  lambda[2] = 1./(hz*hz);

  real_t f = 1.0;
  for (coord_t d = 0; d < 3; d++)
    lambda[d] *= f;

  //matd<real_t> M = Q* (numerics::diag(lambda)*numerics::transpose(Q));
  //symd<real_t> m(3);
  //m.set(M);
  std::pair< vecd<real_t> , matd<real_t> > decomp = {lambda,Q};
  symd<real_t> m(decomp);
  return m;
}

#define MIN(a,b) (a < b) ? (a) : (b);

symd<real_t>
MetricField_UGAWG_Polar2::operator()( const real_t* X ) const {
  real_t eps = 1e-3;
  real_t x[3] = { X[0]+eps, X[1]+eps, X[2] };
  real_t r = std::sqrt( x[0]*x[0] +x[1]*x[1] );
  real_t t = atan2( x[1] , x[0] );

  real_t hz = 0.1;
  real_t ht = 0.1;
  real_t h0 = 1e-3;
  real_t hr = h0 +2.*(0.1 -h0)*fabs( r -0.5 );

  #if 0
  real_t d = 10*(0.6 -r);
  if (d < 0.0) ht = 0.1;
  else ht = d/40. +0.1*(1. -d);
  #else
  real_t d0 = MIN( 10.0 * fabs(r - 0.5) , 1.0 );
  ht = 0.1*d0 + 0.025*(1.0 -d0);
  #endif

  matd<real_t> Q(3,3);
  Q(0,0) = cos(t);
  Q(0,1) = -sin(t);
  Q(0,2) = 0.0;

  Q(1,0) = sin(t);
  Q(1,1) = cos(t);
  Q(1,2) = 0.0;

  Q(2,0) = 0.0;
  Q(2,1) = 0.0;
  Q(2,2) = 1.0;

  vecd<real_t> lambda(3);
  lambda[0] = 1./(hr*hr);
  lambda[1] = 1./(ht*ht);
  lambda[2] = 1./(hz*hz);

  std::pair< vecd<real_t> , matd<real_t> > decomp = {lambda,Q};
  symd<real_t> m(decomp);
  return m;
}

template class MetricField_UniformGeometry<Simplex>;

} // library

} // avro
