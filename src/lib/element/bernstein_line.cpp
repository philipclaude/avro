#include "common/error.h"
#include "element/basis.h"

namespace avro
{

#ifndef ONESIXTH
#define ONESIXTH 0.166666666666666666666666667
#endif

#ifndef ONETHIRD
#define ONETHIRD 0.333333333333333333333333333e+00
#endif

#ifndef ONETWENTYFOURTH
#define ONETWENTYFOURTH 4.1666666666666666666666666666e-02
#endif

/*
 * line, p = 0
*/
template<>
void
Bernstein<Simplex,1,0>::eval( const real_t* x , real_t* phi ) {
  // phi
  phi[0] = 1;
}

template<>
void
Bernstein<Simplex,1,0>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] =  0;
}

template<>
void
Bernstein<Simplex,1,0>::hess( const real_t* x , real_t* phi ) {
  phi[0] =  0;
}

/*
 * line, p = 1
*/
template<>
void
Bernstein<Simplex,1,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] = 1 - s;
  phi[1] =     s;
}

template<>
void
Bernstein<Simplex,1,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = -1;
  phis[1] =  1;
}

template<>
void
Bernstein<Simplex,1,1>::hess( const real_t* x , real_t* phi ) {

  phi[0] = 0;
  phi[1] = 0;
}

/*
 * line, p = 2
*/
template<>
void
Bernstein<Simplex,1,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] = (1 - s)*(1 - s);
  phi[1] = s*s;
  phi[2] = 2*s*(1 - s);
}

template<>
void
Bernstein<Simplex,1,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  real_t* phis = phi;
  phis[0] = 2*(s - 1);
  phis[1] = 2*s;
  phis[2] = 2 - 4*s;
}

template<>
void
Bernstein<Simplex,1,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0] =  2;
  phiss[1] =  2;
  phiss[2] = -4;
}

/*
 * line, p = 3
*/
template<>
void
Bernstein<Simplex,1,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  phi[0] =   t*t*t;
  phi[1] =   s*s*s;
  phi[2] = 3*t*t*s;
  phi[3] = 3*t*s*s;
}

template<>
void
Bernstein<Simplex,1,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  real_t* phis = phi;
  phis[0] = -3*t*t;
  phis[1] =  3*s*s;
  phis[2] =  9*s*s - 12*s + 3;
  phis[3] =  3*(2 - 3*s)*s;
}

template<>
void
Bernstein<Simplex,1,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phiss
  real_t *phiss = phi;
  phiss[0] =  6*t;
  phiss[1] =  6*s;
  phiss[2] = -3*(2 - 3*s);
  phiss[3] =  3*(2 - 3*s);
}

/*
 * line, p = 4
*/
template<>
void
Bernstein<Simplex,1,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phi
  phi[0] =   t*t*t*t;
  phi[1] =   s*s*s*s;
  phi[2] = 6*s*s*t*t;
  phi[3] = 4*s*s*s*t;
  phi[4] = 4*s*t*t*t;
}

template<>
void
Bernstein<Simplex,1,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phis
  real_t* phis = phi;
  phis[0] = -4*t*t*t;
  phis[1] =  4*s*s*s;
  phis[2] =  12*s*(2*s*s - 3*s + 1);
  phis[3] =  4*(3 - 4*s)*s*s;
  phis[4] = -4*t*t*(4*s - 1);
}

template<>
void
Bernstein<Simplex,1,4>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phiss
  real_t* phiss = phi;
  phiss[0] =  12*t*t;
  phiss[1] = 12*s*s;
  phiss[2] = 12*(6*s*s - 6*s + 1);
  phiss[3] = 24*s*(1 - 2*s);
  phiss[4] = -24*(2*s*s - 3*s + 1);
}

/*
 * line, p = 5
*/
template<>
void
Bernstein<Simplex,1,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,1,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,1,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
