//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"
#include "element/basis.h"

namespace avro
{

#ifndef ONESIXTH
#define ONESIXTH 0.16666666666666666666666666666667
#endif

#ifndef ONETHIRD
#define ONETHIRD 0.3333333333333333333333333333333
#endif

#define O3 1./3. // One over 3
#define T3 2./3. // Two over 3

#define Q4 1./4. // Quarter
#define H2 1./2. // Half
#define T4 3./4. // Three quarters

/*
 * triangle, p = 1
*/

template<>
void
Bernstein<Simplex,2,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t b1;
  real_t b2;
  real_t b3;

  b1 = -s + (1 - t);
  b2 = s;
  b3 = 1 - b1 - b2;

  phi[0] = b1;
  phi[1] = b2;
  phi[2] = b3;
}

template<>
void
Bernstein<Simplex,2,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = -1;
  phis[1] = 1;
  phis[2] = 0;

  real_t* phit = phi + 3;
  phit[0] = -1;
  phit[1] = 0;
  phit[2] = 1;
}

template<>
void
Bernstein<Simplex,2,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 9; i++)
    phi[i] = 0;
}

/*
 * triangle, p = 2
*/
template<>
void
Bernstein<Simplex,2,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t b1;
  real_t b2;
  real_t b3;

  b1 = -s + (1 - t);
  b2 = s;
  b3 = 1 - b1 - b2;

  phi[0] = b1*b1;
  phi[1] = b2*b2;
  phi[2] = b3*b3;
  phi[3] = 2*b2*b3;
  phi[4] = 2*b1*b3;
  phi[5] = 2*b1*b2;
}

template<>
void
Bernstein<Simplex,2,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t b1;
  real_t b2;
  real_t b3;

  b1 = -s + (1 - t);
  b2 = s;
  b3 = 1 - b1 - b2;

  real_t* phis = phi;
  phis[0] = -2*b1;
  phis[1] = 2*b2;
  phis[2] = 0;
  phis[3] = 2*b3;
  phis[4] = -2*b3;
  phis[5] = 2*(-b2+b1);

  real_t* phit = phi + 6;
  phit[0] = -2*b1;
  phit[1] = 0;
  phit[2] = 2*b3;
  phit[3] = 2*b2;
  phit[4] = 2*(-b3+b1);
  phit[5] = -2*b2;
}

template<>
void
Bernstein<Simplex,2,2>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangle, p = 3
*/
template<>
void
Bernstein<Simplex,2,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t b1;
  real_t b2;
  real_t b3;

  b1 = -s + (1 - t);
  b2 = s;
  b3 = 1 - b1 - b2;

  phi[0] = b1*b1*b1;
  phi[1] = b2*b2*b2;
  phi[2] = b3*b3*b3;
  phi[3] = 3*b2*b2*b3;
  phi[4] = 3*b2*b3*b3;
  phi[5] = 3*b3*b3*b1;
  phi[6] = 3*b1*b1*b3;
  phi[7] = 3*b1*b1*b2;
  phi[8] = 3*b1*b2*b2;
  phi[9] = 6*b1*b2*b3;
}

template<>
void
Bernstein<Simplex,2,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t* phis = phi;
  phis[0] = -3*(1-s-t)*(1-s-t);
  phis[1] = 3*s*s;
  phis[2] = 0;
  phis[3] = 6*s*t;
  phis[4] = 3*t*t;
  phis[5] = -3*t*t;
  phis[6] = -6*(1-s-t)*t;
  phis[7] = -6*s*(1-s-t)+3*(1-s-t)*(1-s-t);
  phis[8] = -3*s*s+6*s*(1-s-t);
  phis[9] = -6*s*t+6*(1-s-t)*t;

  real_t* phit = phi + 10;
  phit[0] = -3*(1-s-t)*(1-s-t);
  phit[1] = 0;
  phit[2] = 3*t*t;
  phit[3] = 3*s*s;
  phit[4] = 6*s*t;
  phit[5] = 6*(1-s-t)*t-3*t*t;
  phit[6] = 3*(1-s-t)*(1-s-t)-6*(1-s-t)*t;
  phit[7] = -6*s*(1-s-t);
  phit[8] = -3*s*s;
  phit[9] = 6*s*(1-s-t)-6*s*t;
}

template<>
void
Bernstein<Simplex,2,3>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangles, p = 4
*/
template<>
void
Bernstein<Simplex,2,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t b1;
  real_t b2;
  real_t b3;

  b1 = -s + (1 - t);
  b2 = s;
  b3 = 1 - b1 - b2;

  phi[0]  = b1*b1*b1*b1;
  phi[1]  = b2*b2*b2*b2;
  phi[2]  = b3*b3*b3*b3;
  phi[3]  = 4*b2*b2*b2*b3;
  phi[4]  = 6*b2*b2*b3*b3;
  phi[5]  = 4*b3*b3*b3*b2;
  phi[6]  = 4*b3*b3*b3*b1;
  phi[7]  = 6*b3*b3*b1*b1;
  phi[8]  = 4*b1*b1*b1*b3;
  phi[9]  = 4*b1*b1*b1*b2;
  phi[10] = 6*b1*b1*b2*b2;
  phi[11] = 4*b1*b2*b2*b2;
  phi[12] = 12*b1*b1*b2*b3;
  phi[13] = 12*b1*b2*b2*b3;
  phi[14] = 12*b1*b2*b3*b3;
}

template<>
void
Bernstein<Simplex,2,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phis
  real_t* phis = phi;
  phis[0]  = -4*(1-s-t)*(1-s-t)*(1-s-t);
  phis[1]  = 4*s*s*s;
  phis[2]  = 0;
  phis[3]  = 12*s*s*t;
  phis[4]  = 12*s*t*t;
  phis[5]  = 4*t*t*t;
  phis[6]  = -4*t*t*t;
  phis[7]  = -12*(1-s-t)*t*t;
  phis[8]  = -12*(1-s-t)*(1-s-t)*t;
  phis[9]  = -12*s*(1-s-t)*(1-s-t)+4*(1-s-t)*(1-s-t)*(1-s-t) ;
  phis[10] = -12*s*s*(1-s-t)+12*s*(1-s-t)*(1-s-t);
  phis[11] = -4*s*s*s+12*s*s*(1-s-t);
  phis[12] = -24*s*(1-s-t)*t+12*(1-s-t)*(1-s-t)*t;
  phis[13] = -12*s*s*t+24*s*(1-s-t)*t;
  phis[14] = -12*s*t*t+12*(1-s-t)*t*t;

  // phit
  real_t *phit = phi + 15;
  phit[0]  = -4*(1-s-t)*(1-s-t)*(1-s-t);
  phit[1]  = 0;
  phit[2]  = 4*t*t*t;
  phit[3]  = 4*s*s*s;
  phit[4]  = 12*s*s*t;
  phit[5]  = 12*s*t*t;
  phit[6]  = 12*(1-s-t)*t*t-4*t*t*t;
  phit[7]  = 12*(1-s-t)*(1-s-t)*t-12*(1-s-t)*t*t;
  phit[8]  = 4*(1-s-t)*(1-s-t)*(1-s-t)-12*(1-s-t)*(1-s-t)*t;
  phit[9]  = -12*s*(1-s-t)*(1-s-t);
  phit[10] = -12*s*s*(1-s-t);
  phit[11] = -4*s*s*s;
  phit[12] = 12*s*(1-s-t)*(1-s-t)-24*s*(1-s-t)*t;
  phit[13] = 12*s*s*(1-s-t)-12*s*s*t;
  phit[14] = 24*s*(1-s-t)*t-12*s*t*t;
}

template<>
void
Bernstein<Simplex,2,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangles, p = 5
*/
template<>
void
Bernstein<Simplex,2,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,2,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,2,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
