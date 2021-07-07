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
Lagrange<Simplex,2,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phi
  phi[0] =  -s-t+1.0;
  phi[1] =  s;
  phi[2] =  t;
}

template<>
void
Lagrange<Simplex,2,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] =  -1.0;
  phis[1] =  1.0;
  phis[2] =  0.0;

  real_t* phit = phi + 3;
  phit[0] =  -1.0;
  phit[1] =  0.0;
  phit[2] =  1.0;
}

template<>
void
Lagrange<Simplex,2,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 9; i++)
    phi[i] = 0;
}

/*
 * triangle, p = 2
*/
template<>
void
Lagrange<Simplex,2,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  phi[0] =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
  phi[1] =  -s+(s*s)*2.0;
  phi[2] =  -t+(t*t)*2.0;
  phi[3] =  s*t*4.0;
  phi[4] =  t*(s+t-1.0)*-4.0;
  phi[5] =  -s*(s*4.0+t*4.0-4.0);
}

template<>
void
Lagrange<Simplex,2,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t* phis = phi;
  phis[0] =  s*4.0+t*4.0-3.0;
  phis[1] =  s*4.0-1.0;
  phis[2] =  0.0;
  phis[3] =  t*4.0;
  phis[4] =  t*-4.0;
  phis[5] =  s*-8.0-t*4.0+4.0;

  real_t* phit = phi + 6;
  phit[0] =  s*4.0+t*4.0-3.0;
  phit[1] =  0.0;
  phit[2] =  t*4.0-1.0;
  phit[3] =  s*4.0;
  phit[4] =  s*-4.0-t*8.0+4.0;
  phit[5] =  s*-4.0;
}

template<>
void
Lagrange<Simplex,2,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0] =  4.0;
  phiss[1] =  4.0;
  phiss[2] =  0.0;
  phiss[3] =  0.0;
  phiss[4] =  0.0;
  phiss[5] =  -8.0;

  // phist
  real_t *phist = phi + 6;
  phist[0] =  4.0;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  4.0;
  phist[4] =  -4.0;
  phist[5] =  -4.0;

  // phitt
  real_t *phitt = phi + 12;
  phitt[0] =  4.0;
  phitt[1] =  0.0;
  phitt[2] =  4.0;
  phitt[3] =  0.0;
  phitt[4] =  -8.0;
  phitt[5] =  0.0;
}

/*
 * triangle, p = 3
*/
template<>
void
Lagrange<Simplex,2,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
              (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;
  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);
  phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);
  phi[3] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);
  phi[4] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);
  phi[5] =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);
  phi[6] =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);
  phi[7] =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);
  phi[8] =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);
  phi[9] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;
}

template<>
void
Lagrange<Simplex,2,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t* phis = phi;
  phis[0] =  s*1.8E1+t*1.8E1-s*t*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-1.1E1/2.0;
  phis[1] =  s*-9.0+(s*s)*(2.7E1/2.0)+1.0;
  phis[2] =  0.0;
  phis[3] =  t*(-9.0/2.0)+s*t*2.7E1;
  phis[4] =  t*(-9.0/2.0)+(t*t)*(2.7E1/2.0);
  phis[5] =  t*(9.0/2.0)-(t*t)*(2.7E1/2.0);
  phis[6] =  t*(-4.5E1/2.0)+s*t*2.7E1+(t*t)*2.7E1;
  phis[7] =  s*-4.5E1-t*(4.5E1/2.0)+s*t*5.4E1+(s*s)*(8.1E1/2.0)+(t*t)*(2.7E1/2.0)+9.0;
  phis[8] =  s*3.6E1+t*(9.0/2.0)-s*t*2.7E1-(s*s)*(8.1E1/2.0)-9.0/2.0;
  phis[9] =  t*2.7E1-s*t*5.4E1-(t*t)*2.7E1;

  real_t* phit = phi + 10;
  phit[0] =  s*1.8E1+t*1.8E1-s*t*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-1.1E1/2.0;
  phit[1] =  0.0;
  phit[2] =  t*-9.0+(t*t)*(2.7E1/2.0)+1.0;
  phit[3] =  s*(-9.0/2.0)+(s*s)*(2.7E1/2.0);
  phit[4] =  s*(-9.0/2.0)+s*t*2.7E1;
  phit[5] =  s*(9.0/2.0)+t*3.6E1-s*t*2.7E1-(t*t)*(8.1E1/2.0)-9.0/2.0;
  phit[6] =  s*(-4.5E1/2.0)-t*4.5E1+s*t*5.4E1+(s*s)*(2.7E1/2.0)+(t*t)*(8.1E1/2.0)+9.0;
  phit[7] =  s*(-4.5E1/2.0)+s*t*2.7E1+(s*s)*2.7E1;
  phit[8] =  s*(9.0/2.0)-(s*s)*(2.7E1/2.0);
  phit[9] =  s*2.7E1-s*t*5.4E1-(s*s)*2.7E1;
}

template<>
void
Lagrange<Simplex,2,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phiss
  real_t *phiss = phi;
  phiss[0] =  s*-2.7E1-t*2.7E1+1.8E1;
  phiss[1] =  s*2.7E1-9.0;
  phiss[2] =  0.0;
  phiss[3] =  t*2.7E1;
  phiss[4] =  0.0;
  phiss[5] =  0.0;
  phiss[6] =  t*2.7E1;
  phiss[7] =  s*8.1E1+t*5.4E1-4.5E1;
  phiss[8] =  s*-8.1E1-t*2.7E1+3.6E1;
  phiss[9] =  t*-5.4E1;

  // phist
  real_t *phist = phi + 10;
  phist[0] =  s*-2.7E1-t*2.7E1+1.8E1;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  s*2.7E1-9.0/2.0;
  phist[4] =  t*2.7E1-9.0/2.0;
  phist[5] =  t*-2.7E1+9.0/2.0;
  phist[6] =  s*2.7E1+t*5.4E1-4.5E1/2.0;
  phist[7] =  s*5.4E1+t*2.7E1-4.5E1/2.0;
  phist[8] =  s*-2.7E1+9.0/2.0;
  phist[9] =  s*-5.4E1-t*5.4E1+2.7E1;

  // phitt
  real_t* phitt = phi + 20;
  phitt[0] =  s*-2.7E1-t*2.7E1+1.8E1;
  phitt[1] =  0.0;
  phitt[2] =  t*2.7E1-9.0;
  phitt[3] =  0.0;
  phitt[4] =  s*2.7E1;
  phitt[5] =  s*-2.7E1-t*8.1E1+3.6E1;
  phitt[6] =  s*5.4E1+t*8.1E1-4.5E1;
  phitt[7] =  s*2.7E1;
  phitt[8] =  0.0;
  phitt[9] =  s*-5.4E1;
}

/*
 * triangles, p = 4
*/
template<>
void
Lagrange<Simplex,2,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phi
  phi[0] =  (s*s)*(t*t)*6.4E1-ONETHIRD*s*2.5E1-ONETHIRD*t*2.5E1+ONETHIRD*(s*s)*7.0E1-ONETHIRD*(s*s*s)*8.0E1+
      ONETHIRD*(s*s*s*s)*3.2E1+ONETHIRD*(t*t)*7.0E1-ONETHIRD*(t*t*t)*8.0E1+ONETHIRD*(t*t*t*t)*3.2E1-s*(t*t)*
      8.0E1-(s*s)*t*8.0E1+ONETHIRD*s*t*1.4E2+ONETHIRD*s*(t*t*t)*1.28E2+ONETHIRD*(s*s*s)*t*1.28E2+1.0;
  phi[1] =  -s+ONETHIRD*(s*s)*2.2E1+ONETHIRD*(s*s*s*s)*3.2E1-(s*s*s)*1.6E1;
  phi[2] =  -t+ONETHIRD*(t*t)*2.2E1+ONETHIRD*(t*t*t*t)*3.2E1-(t*t*t)*1.6E1;
  phi[3] =  (s*s)*t*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*(s*s*s)*t*1.28E2;
  phi[4] =  (s*s)*(t*t)*6.4E1+s*t*4.0-s*(t*t)*1.6E1-(s*s)*t*1.6E1;
  phi[5] =  s*(t*t)*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*s*(t*t*t)*1.28E2;
  phi[6] =  ONETHIRD*t*1.6E1-ONETHIRD*(t*t)*1.12E2+ONETHIRD*(t*t*t)*2.24E2-ONETHIRD*(t*t*t*t)*1.28E2+s*
      (t*t)*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*s*(t*t*t)*1.28E2;
  phi[7] =  t*-1.2E1+(s*s)*(t*t)*6.4E1+s*t*2.8E1-s*(t*t)*1.44E2-(s*s)*t*1.6E1+s*(t*t*t)*1.28E2+(t*t)*7.6E1-
      (t*t*t)*1.28E2+(t*t*t*t)*6.4E1;
  phi[8] =  t*1.6E1-(s*s)*(t*t)*1.28E2-ONETHIRD*(t*t)*2.08E2-ONETHIRD*(t*t*t*t)*1.28E2+s*(t*t)*1.92E2+(s*
      s)*t*9.6E1-s*(t*t*t)*1.28E2+(t*t*t)*9.6E1-ONETHIRD*s*t*2.08E2-ONETHIRD*(s*s*s)*t*1.28E2;
  phi[9] =  s*1.6E1-(s*s)*(t*t)*1.28E2-ONETHIRD*(s*s)*2.08E2-ONETHIRD*(s*s*s*s)*1.28E2+s*(t*t)*9.6E1+(s*
      s)*t*1.92E2-(s*s*s)*t*1.28E2+(s*s*s)*9.6E1-ONETHIRD*s*t*2.08E2-ONETHIRD*s*(t*t*t)*1.28E2;
  phi[10] =  s*-1.2E1+(s*s)*(t*t)*6.4E1+s*t*2.8E1-s*(t*t)*1.6E1-(s*s)*t*1.44E2+(s*s*s)*t*1.28E2+(s*s)*7.6E1-
      (s*s*s)*1.28E2+(s*s*s*s)*6.4E1;
  phi[11] =  ONETHIRD*s*1.6E1-ONETHIRD*(s*s)*1.12E2+ONETHIRD*(s*s*s)*2.24E2-ONETHIRD*(s*s*s*s)*1.28E2+(s*
      s)*t*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*(s*s*s)*t*1.28E2;
  phi[12] =  (s*s)*(t*t)*2.56E2+s*t*9.6E1-s*(t*t)*2.24E2-(s*s)*t*2.24E2+s*(t*t*t)*1.28E2+(s*s*s)*t*1.28E2;
  phi[13] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*3.2E1+(s*s)*t*1.6E2-(s*s*s)*t*1.28E2;
  phi[14] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*1.6E2+(s*s)*t*3.2E1-s*(t*t*t)*1.28E2;
}

template<>
void
Lagrange<Simplex,2,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phis
  real_t* phis = phi;
  phis[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2-s*t*1.6E2-ONETHIRD*(s*s)*2.4E2+ONETHIRD*
      (s*s*s)*1.28E2+ONETHIRD*(t*t*t)*1.28E2+s*(t*t)*1.28E2-(t*t)*8.0E1+ONETHIRD*(s*s)*t*3.84E2;
  phis[1] =  ONETHIRD*s*4.4E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*4.8E1-1.0;
  phis[2] =  0.0;
  phis[3] =  ONETHIRD*t*1.6E1-s*t*6.4E1+ONETHIRD*(s*s)*t*3.84E2;
  phis[4] =  t*4.0-s*t*3.2E1+s*(t*t)*1.28E2-(t*t)*1.6E1;
  phis[5] =  ONETHIRD*t*1.6E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*3.2E1;
  phis[6] =  ONETHIRD*t*-1.6E1-ONETHIRD*(t*t*t)*1.28E2+(t*t)*3.2E1;
  phis[7] =  t*2.8E1-s*t*3.2E1+s*(t*t)*1.28E2-(t*t)*1.44E2+(t*t*t)*1.28E2;
  phis[8] =  ONETHIRD*t*-2.08E2+s*t*1.92E2-s*(t*t)*2.56E2+(t*t)*1.92E2-(t*t*t)*1.28E2-ONETHIRD*(s*s)*t*
      3.84E2;
  phis[9] =  ONETHIRD*s*-4.16E2-ONETHIRD*t*2.08E2+s*t*3.84E2-ONETHIRD*(s*s*s)*5.12E2-ONETHIRD*(t*t*t)*1.28E2-
      s*(t*t)*2.56E2-(s*s)*t*3.84E2+(s*s)*2.88E2+(t*t)*9.6E1+1.6E1;
  phis[10] =  s*1.52E2+t*2.8E1-s*t*2.88E2+s*(t*t)*1.28E2+(s*s)*t*3.84E2-(s*s)*3.84E2+(s*s*s)*2.56E2-(t*
      t)*1.6E1-1.2E1;
  phis[11] =  ONETHIRD*1.6E1-ONETHIRD*s*2.24E2-ONETHIRD*t*1.6E1+s*t*6.4E1+ONETHIRD*(s*s)*6.72E2-ONETHIRD*
      (s*s*s)*5.12E2-ONETHIRD*(s*s)*t*3.84E2;
  phis[12] =  t*9.6E1-s*t*4.48E2+s*(t*t)*5.12E2+(s*s)*t*3.84E2-(t*t)*2.24E2+(t*t*t)*1.28E2;
  phis[13] =  t*-3.2E1+s*t*3.2E2-s*(t*t)*2.56E2-(s*s)*t*3.84E2+(t*t)*3.2E1;
  phis[14] =  t*-3.2E1+s*t*6.4E1-s*(t*t)*2.56E2+(t*t)*1.6E2-(t*t*t)*1.28E2;

  // phit
  real_t *phit = phi + 15;
  phit[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2-s*t*1.6E2+ONETHIRD*(s*s*s)*1.28E2-ONETHIRD*
      (t*t)*2.4E2+ONETHIRD*(t*t*t)*1.28E2+(s*s)*t*1.28E2-(s*s)*8.0E1+ONETHIRD*s*(t*t)*3.84E2;
  phit[1] =  0.0;
  phit[2] =  ONETHIRD*t*4.4E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*4.8E1-1.0;
  phit[3] =  ONETHIRD*s*1.6E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*3.2E1;
  phit[4] =  s*4.0-s*t*3.2E1+(s*s)*t*1.28E2-(s*s)*1.6E1;
  phit[5] =  ONETHIRD*s*1.6E1-s*t*6.4E1+ONETHIRD*s*(t*t)*3.84E2;
  phit[6] =  ONETHIRD*1.6E1-ONETHIRD*s*1.6E1-ONETHIRD*t*2.24E2+s*t*6.4E1+ONETHIRD*(t*t)*6.72E2-ONETHIRD*
      (t*t*t)*5.12E2-ONETHIRD*s*(t*t)*3.84E2;
  phit[7] =  s*2.8E1+t*1.52E2-s*t*2.88E2+s*(t*t)*3.84E2+(s*s)*t*1.28E2-(s*s)*1.6E1-(t*t)*3.84E2+(t*t*t)*
      2.56E2-1.2E1;
  phit[8] =  ONETHIRD*s*-2.08E2-ONETHIRD*t*4.16E2+s*t*3.84E2-ONETHIRD*(s*s*s)*1.28E2-ONETHIRD*(t*t*t)*5.12E2-
      s*(t*t)*3.84E2-(s*s)*t*2.56E2+(s*s)*9.6E1+(t*t)*2.88E2+1.6E1;
  phit[9] =  ONETHIRD*s*-2.08E2+s*t*1.92E2-(s*s)*t*2.56E2+(s*s)*1.92E2-(s*s*s)*1.28E2-ONETHIRD*s*(t*t)*
      3.84E2;
  phit[10] =  s*2.8E1-s*t*3.2E1+(s*s)*t*1.28E2-(s*s)*1.44E2+(s*s*s)*1.28E2;
  phit[11] =  ONETHIRD*s*-1.6E1-ONETHIRD*(s*s*s)*1.28E2+(s*s)*3.2E1;
  phit[12] =  s*9.6E1-s*t*4.48E2+s*(t*t)*3.84E2+(s*s)*t*5.12E2-(s*s)*2.24E2+(s*s*s)*1.28E2;
  phit[13] =  s*-3.2E1+s*t*6.4E1-(s*s)*t*2.56E2+(s*s)*1.6E2-(s*s*s)*1.28E2;
  phit[14] =  s*-3.2E1+s*t*3.2E2-s*(t*t)*3.84E2-(s*s)*t*2.56E2+(s*s)*3.2E1;
}

template<>
void
Lagrange<Simplex,2,4>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phiss
  real_t* phiss = phi;
  phiss[0] =  ONETHIRD*1.4E2-t*1.6E2-ONETHIRD*s*4.8E2+ONETHIRD*(s*s)*3.84E2+(t*t)*1.28E2+ONETHIRD*s*t*7.68E2;
  phiss[1] =  ONETHIRD*4.4E1-s*9.6E1+ONETHIRD*(s*s)*3.84E2;
  phiss[2] =  0.0;
  phiss[3] =  t*-6.4E1+ONETHIRD*s*t*7.68E2;
  phiss[4] =  t*-3.2E1+(t*t)*1.28E2;
  phiss[5] =  0.0;
  phiss[6] =  0.0;
  phiss[7] =  t*-3.2E1+(t*t)*1.28E2;
  phiss[8] =  t*1.92E2-(t*t)*2.56E2-ONETHIRD*s*t*7.68E2;
  phiss[9] =  ONETHIRD*-4.16E2+s*5.76E2+t*3.84E2-s*t*7.68E2-ONETHIRD*(s*s)*1.536E3-(t*t)*2.56E2;
  phiss[10] =  s*-7.68E2-t*2.88E2+s*t*7.68E2+(s*s)*7.68E2+(t*t)*1.28E2+1.52E2;
  phiss[11] =  ONETHIRD*-2.24E2+t*6.4E1+ONETHIRD*s*1.344E3-ONETHIRD*(s*s)*1.536E3-ONETHIRD*s*t*7.68E2;
  phiss[12] =  t*-4.48E2+s*t*7.68E2+(t*t)*5.12E2;
  phiss[13] =  t*3.2E2-s*t*7.68E2-(t*t)*2.56E2;
  phiss[14] =  t*6.4E1-(t*t)*2.56E2;

  // phist
  real_t* phist = phi + 15;
  phist[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2+s*t*2.56E2+ONETHIRD*(s*s)*3.84E2+ONETHIRD*(t*t)*3.84E2;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  ONETHIRD*1.6E1-s*6.4E1+ONETHIRD*(s*s)*3.84E2;
  phist[4] =  s*-3.2E1-t*3.2E1+s*t*2.56E2+4.0;
  phist[5] =  ONETHIRD*1.6E1-t*6.4E1+ONETHIRD*(t*t)*3.84E2;
  phist[6] =  ONETHIRD*-1.6E1+t*6.4E1-ONETHIRD*(t*t)*3.84E2;
  phist[7] =  s*-3.2E1-t*2.88E2+s*t*2.56E2+(t*t)*3.84E2+2.8E1;
  phist[8] =  ONETHIRD*-2.08E2+s*1.92E2+t*3.84E2-s*t*5.12E2-ONETHIRD*(s*s)*3.84E2-(t*t)*3.84E2;
  phist[9] =  ONETHIRD*-2.08E2+s*3.84E2+t*1.92E2-s*t*5.12E2-ONETHIRD*(t*t)*3.84E2-(s*s)*3.84E2;
  phist[10] =  s*-2.88E2-t*3.2E1+s*t*2.56E2+(s*s)*3.84E2+2.8E1;
  phist[11] =  ONETHIRD*-1.6E1+s*6.4E1-ONETHIRD*(s*s)*3.84E2;
  phist[12] =  s*-4.48E2-t*4.48E2+s*t*1.024E3+(s*s)*3.84E2+(t*t)*3.84E2+9.6E1;
  phist[13] =  s*3.2E2+t*6.4E1-s*t*5.12E2-(s*s)*3.84E2-3.2E1;
  phist[14] =  s*6.4E1+t*3.2E2-s*t*5.12E2-(t*t)*3.84E2-3.2E1;

  // phitt
  real_t* phitt = phi + 30;
  phitt[0] =  ONETHIRD*1.4E2-s*1.6E2-ONETHIRD*t*4.8E2+ONETHIRD*(t*t)*3.84E2+(s*s)*1.28E2+ONETHIRD*s*t*7.68E2;
  phitt[1] =  0.0;
  phitt[2] =  ONETHIRD*4.4E1-t*9.6E1+ONETHIRD*(t*t)*3.84E2;
  phitt[3] =  0.0;
  phitt[4] =  s*-3.2E1+(s*s)*1.28E2;
  phitt[5] =  s*-6.4E1+ONETHIRD*s*t*7.68E2;
  phitt[6] =  ONETHIRD*-2.24E2+s*6.4E1+ONETHIRD*t*1.344E3-ONETHIRD*(t*t)*1.536E3-ONETHIRD*s*t*7.68E2;
  phitt[7] =  s*-2.88E2-t*7.68E2+s*t*7.68E2+(s*s)*1.28E2+(t*t)*7.68E2+1.52E2;
  phitt[8] =  ONETHIRD*-4.16E2+s*3.84E2+t*5.76E2-s*t*7.68E2-ONETHIRD*(t*t)*1.536E3-(s*s)*2.56E2;
  phitt[9] =  s*1.92E2-(s*s)*2.56E2-ONETHIRD*s*t*7.68E2;
  phitt[10] =  s*-3.2E1+(s*s)*1.28E2;
  phitt[11] =  0.0;
  phitt[12] =  s*-4.48E2+s*t*7.68E2+(s*s)*5.12E2;
  phitt[13] =  s*6.4E1-(s*s)*2.56E2;
  phitt[14] =  s*3.2E2-s*t*7.68E2-(s*s)*2.56E2;
}

/*
 * triangles, p = 5
*/
template<>
void
Lagrange<Simplex,2,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,2,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,2,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
