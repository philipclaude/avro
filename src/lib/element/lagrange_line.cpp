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
 * line, p = 1
*/

template<>
void
Lagrange<Simplex,1,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] =  -s+1.0;
  phi[1] =  s;
}

template<>
void
Lagrange<Simplex,1,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] =  -1.0;
  phis[1] =  1.0;
}

template<>
void
Lagrange<Simplex,1,1>::hess( const real_t* x , real_t* phi ) {

  phi[0] =  0.0;
  phi[1] =  0.0;
}

/*
 * line, p = 2
*/
template<>
void
Lagrange<Simplex,1,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] =  s*-3.0+(s*s)*2.0+1.0;
  phi[1] =  -s+(s*s)*2.0;
  phi[2] =  s*4.0-(s*s)*4.0;
}

template<>
void
Lagrange<Simplex,1,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  real_t* phis = phi;
  phis[0] =  s*4.0-3.0;
  phis[1] =  s*4.0-1.0;
  phis[2] =  s*-8.0+4.0;
}

template<>
void
Lagrange<Simplex,1,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0] =  4.0;
  phiss[1] =  4.0;
  phiss[2] =  -8.0;
}

/*
 * line, p = 3
*/
template<>
void
Lagrange<Simplex,1,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  phi[0] =  s*(-1.1E1/2.0)+(s*s)*9.0-(s*s*s)*(9.0/2.0)+1.0;
  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);
  phi[2] =  s*9.0-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);
  phi[3] =  s*(-9.0/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);
}

template<>
void
Lagrange<Simplex,1,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  real_t* phis = phi;
  phis[0] =  s*1.8E1-(s*s)*(2.7E1/2.0)-1.1E1/2.0;
  phis[1] =  s*-9.0+(s*s)*(2.7E1/2.0)+1.0;
  phis[2] =  s*-4.5E1+(s*s)*(8.1E1/2.0)+9.0;
  phis[3] =  s*3.6E1-(s*s)*(8.1E1/2.0)-9.0/2.0;
}

template<>
void
Lagrange<Simplex,1,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phiss
  real_t *phiss = phi;
  phiss[0] =  s*-2.7E1+1.8E1;
  phiss[1] =  s*2.7E1-9.0;
  phiss[2] =  s*8.1E1-4.5E1;
  phiss[3] =  s*-8.1E1+3.6E1;
}

/*
 * line, p = 4
*/
template<>
void
Lagrange<Simplex,1,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phi
  phi[0] =  ONETHIRD*(s*-2.5E1+(s*s)*7.0E1-(s*s*s)*8.0E1+(s*s*s*s)*3.2E1+3.0);
  phi[1] =  -ONETHIRD*(s*3.0-(s*s)*2.2E1+(s*s*s)*4.8E1-(s*s*s*s)*3.2E1);
  phi[2] =  ONETHIRD*(s*4.8E1-(s*s)*2.08E2+(s*s*s)*2.88E2-(s*s*s*s)*1.28E2);
  phi[3] =  -ONETHIRD*(s*3.6E1-(s*s)*2.28E2+(s*s*s)*3.84E2-(s*s*s*s)*1.92E2);
  phi[4] =  ONETHIRD*(s*1.6E1-(s*s)*1.12E2+(s*s*s)*2.24E2-(s*s*s*s)*1.28E2);
}

template<>
void
Lagrange<Simplex,1,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phis
  real_t* phis = phi;
  phis[0] =  ONETHIRD*(s*1.4E2-(s*s)*2.4E2+(s*s*s)*1.28E2-2.5E1);
  phis[1] =  ONETHIRD*(s*4.4E1-(s*s)*1.44E2+(s*s*s)*1.28E2-3.0);
  phis[2] =  -ONETHIRD*(s*4.16E2-(s*s)*8.64E2+(s*s*s)*5.12E2-4.8E1);
  phis[3] =  ONETHIRD*(s*4.56E2-(s*s)*1.152E3+(s*s*s)*7.68E2-3.6E1);
  phis[4] =  -ONETHIRD*(s*2.24E2-(s*s)*6.72E2+(s*s*s)*5.12E2-1.6E1);
}

template<>
void
Lagrange<Simplex,1,4>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phiss
  real_t* phiss = phi;
  phiss[0] =  ONETHIRD*(s*-4.8E2+(s*s)*3.84E2+1.4E2);
  phiss[1] =  ONETHIRD*(s*-2.88E2+(s*s)*3.84E2+4.4E1);
  phiss[2] =  -ONETHIRD*(s*-1.728E3+(s*s)*1.536E3+4.16E2);
  phiss[3] =  ONETHIRD*(s*-2.304E3+(s*s)*2.304E3+4.56E2);
  phiss[4] =  -ONETHIRD*(s*-1.344E3+(s*s)*1.536E3+2.24E2);
}

/*
 * line, p = 5
*/
template<>
void
Lagrange<Simplex,1,5>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  phi[0] =  -ONETWENTYFOURTH*(s*2.74E2-(s*s)*1.125E3+(s*s*s)*2.125E3-(s*s*s*s)*1.875E3+(s*s*s*s*s)*6.25E2-2.4E1);
  phi[1] =  ONETWENTYFOURTH*(s*2.4E1-(s*s)*2.5E2+(s*s*s)*8.75E2-(s*s*s*s)*1.25E3+(s*s*s*s*s)*6.25E2);
  phi[2] =  ONETWENTYFOURTH*(s*6.0E2-(s*s)*3.85E3+(s*s*s)*8.875E3-(s*s*s*s)*8.75E3+(s*s*s*s*s)*3.125E3);
  phi[3] =  -ONETWENTYFOURTH*(s*6.0E2-(s*s)*5.35E3+(s*s*s)*1.475E4-(s*s*s*s)*1.625E4+(s*s*s*s*s)*6.25E3);
  phi[4] =  ONETWENTYFOURTH*(s*4.0E2-(s*s)*3.9E3+(s*s*s)*1.225E4-(s*s*s*s)*1.5E4+(s*s*s*s*s)*6.25E3);
  phi[5] =  -ONETWENTYFOURTH*(s*1.5E2-(s*s)*1.525E3+(s*s*s)*5.125E3-(s*s*s*s)*6.875E3+(s*s*s*s*s)*3.125E3);
}

template<>
void
Lagrange<Simplex,1,5>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phis
  real_t* phis = phi;
  phis[0] =  -ONETWENTYFOURTH*(s*-2.25E3+(s*s)*6.375E3-(s*s*s)*7.5E3+(s*s*s*s)*3.125E3+2.74E2);
  phis[1] =  ONETWENTYFOURTH*(s*-5.0E2+(s*s)*2.625E3-(s*s*s)*5.0E3+(s*s*s*s)*3.125E3+2.4E1);
  phis[2] =  ONETWENTYFOURTH*(s*-7.7E3+(s*s)*2.6625E4-(s*s*s)*3.5E4+(s*s*s*s)*1.5625E4+6.0E2);
  phis[3] =  -ONETWENTYFOURTH*(s*-1.07E4+(s*s)*4.425E4-(s*s*s)*6.5E4+(s*s*s*s)*3.125E4+6.0E2);
  phis[4] =  ONETWENTYFOURTH*(s*-7.8E3+(s*s)*3.675E4-(s*s*s)*6.0E4+(s*s*s*s)*3.125E4+4.0E2);
  phis[5] =  -ONETWENTYFOURTH*(s*-3.05E3+(s*s)*1.5375E4-(s*s*s)*2.75E4+(s*s*s*s)*1.5625E4+1.5E2);
}

template<>
void
Lagrange<Simplex,1,5>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];

  // phiss
  real_t* phiss = phi;
  phiss[0] =  -ONETWENTYFOURTH*(s*1.275E4-(s*s)*2.25E4+(s*s*s)*1.25E4-2.25E3);
  phiss[1] =  ONETWENTYFOURTH*(s*5.25E3-(s*s)*1.5E4+(s*s*s)*1.25E4-5.0E2);
  phiss[2] =  ONETWENTYFOURTH*(s*5.325E4-(s*s)*1.05E5+(s*s*s)*6.25E4-7.7E3);
  phiss[3] =  -ONETWENTYFOURTH*(s*8.85E4-(s*s)*1.95E5+(s*s*s)*1.25E5-1.07E4);
  phiss[4] =  ONETWENTYFOURTH*(s*7.35E4-(s*s)*1.8E5+(s*s*s)*1.25E5-7.8E3);
  phiss[5] =  -ONETWENTYFOURTH*(s*3.075E4-(s*s)*8.25E4+(s*s*s)*6.25E4-3.05E3);
}

} // avro
