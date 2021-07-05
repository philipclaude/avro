#include "common/error.h"
#include "element/basis.h"

#include <cmath>

namespace avro
{

/*
 * line, p = 0
*/
template<>
void
Legendre<Simplex,1,0>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] =  1;
}

template<>
void
Legendre<Simplex,1,0>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = 0;
}

template<>
void
Legendre<Simplex,1,0>::hess( const real_t* x , real_t* phi ) {
  phi[0] = 0;
}

/*
 * line, p = 1
*/
template<>
void
Legendre<Simplex,1,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  real_t t = 1 - s;

  phi[0] = 1;
  phi[1] = sqrt(3)*(s - t);
}

template<>
void
Legendre<Simplex,1,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] =  0;
  phis[1] =  2*sqrt(3);
}

template<>
void
Legendre<Simplex,1,1>::hess( const real_t* x , real_t* phi ) {

  phi[0] =  0.0;
  phi[1] =  0.0;
}

/*
 * line, p = 2
*/
template<>
void
Legendre<Simplex,1,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  real_t t = 1 - s;

  phi[0] = 1;
  phi[1] = sqrt(3)*(s - t);
  phi[2] = sqrt(5)*(1 - 6*s*t);
}

template<>
void
Legendre<Simplex,1,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  real_t* phis = phi;
  phis[0] =  0;
  phis[1] =  2*sqrt(3);
  phis[2] =  6*sqrt(5)*(s - t);
}

template<>
void
Legendre<Simplex,1,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0] = 0;
  phiss[1] = 0;
  phiss[2] = 12*sqrt(5);
}

/*
 * line, p = 3
*/
template<>
void
Legendre<Simplex,1,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  phi[0] = 1;
  phi[1] = sqrt(3)*(s - t);
  phi[2] = sqrt(5)*(1 - 6*s*t);
  phi[3] = sqrt(7)*(1 - 10*s*t)*(s - t);
}

template<>
void
Legendre<Simplex,1,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  real_t* phis = phi;
  phis[0] =  0;
  phis[1] =  2*sqrt(3);
  phis[2] =  6*sqrt(5)*(s - t);
  phis[3] = 12*sqrt(7)*(1 - 5*s*t);
}

template<>
void
Legendre<Simplex,1,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phiss
  real_t *phiss = phi;

  phiss[0] = 0;
  phiss[1] = 0;
  phiss[2] = 12*sqrt(5);
  phiss[3] = 60*sqrt(7)*(s - t);
}

/*
 * line, p = 4
*/
template<>
void
Legendre<Simplex,1,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phi
  phi[0] = 1;
  phi[1] = sqrt(3)*(s - t);
  phi[2] = sqrt(5)*(1 - 6*s*t);
  phi[3] = sqrt(7)*(1 - 10*s*t)*(s - t);
  phi[4] = 3*(1 - 10*s*t*(2 - 7*s*t));
}

template<>
void
Legendre<Simplex,1,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phis
  real_t* phis = phi;
  phis[0] =  0;
  phis[1] =  2*sqrt(3);
  phis[2] =  6*sqrt(5)*(s - t);
  phis[3] = 12*sqrt(7)*(1 - 5*s*t);
  phis[4] = 60*(1 - 7*s*t)*(s - t);
}

template<>
void
Legendre<Simplex,1,4>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phiss
  real_t* phiss = phi;
  phiss[0] = 0;
  phiss[1] = 0;
  phiss[2] = 12*sqrt(5);
  phiss[3] = 60*sqrt(7)*(s - t);
  phiss[4] = 180*(3 - 14*s*t);
}

/*
 * line, p = 5
*/
template<>
void
Legendre<Simplex,1,5>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  phi[0] = 1;
  phi[1] = sqrt(3)*(s - t);
  phi[2] = sqrt(5)*(1 - 6*s*t);
  phi[3] = sqrt(7)*(1 - 10*s*t)*(s - t);
  phi[4] = 3*(1 - 10*s*t*(2 - 7*s*t));
  phi[5] = sqrt(11)*(1 - 14*s*t*(s - 2*t)*(2*s - t))*(s - t);
}

template<>
void
Legendre<Simplex,1,5>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phis
  real_t* phis = phi;
  phis[0] =  0;
  phis[1] =  2*sqrt(3);
  phis[2] =  6*sqrt(5)*(s - t);
  phis[3] = 12*sqrt(7)*(1 - 5*s*t);
  phis[4] = 60*(1 - 7*s*t)*(s - t);
  phis[5] = 30*sqrt(11)*(1 - 14*s*t*(1 - 3*s*t));
}

template<>
void
Legendre<Simplex,1,5>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = 1 - s;

  // phiss
  real_t* phiss = phi;
  phiss[0] = 0;
  phiss[1] = 0;
  phiss[2] = 12*sqrt(5);
  phiss[3] = 60*sqrt(7)*(s - t);
  phiss[4] = 180*(3 - 14*s*t);
  phiss[5] = 420*sqrt(11)*(1 - 6*s*t)*(s - t);
}

} // avro
