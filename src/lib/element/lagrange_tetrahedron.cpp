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
 * pentatopes, p = 1
*/

template<>
void
Lagrange<Simplex,3,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  phi[0] = 1 - s - t - u; // 1 @ node 0 (s = 0, t = 0, u = 0)
  phi[1] =     s        ; // 1 @ node 1 (s = 1, t = 0, u = 0)
  phi[2] =         t    ; // 1 @ node 2 (s = 0, t = 1, u = 0)
  phi[3] =             u; // 1 @ node 3 (s = 0, t = 0, u = 1)
}

template<>
void
Lagrange<Simplex,3,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = -1;
  phis[1] =  1;
  phis[2] =  0;
  phis[3] =  0;

  real_t* phit = phi + 4;
  phit[0] = -1;
  phit[1] =  0;
  phit[2] =  1;
  phit[3] =  0;

  real_t* phiu = phi + 8;
  phiu[0]= -1;
  phiu[1]= 0;
  phiu[2]= 0;
  phiu[3]= 1;
}

template<>
void
Lagrange<Simplex,3,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 24; i++)
    phi[i] = 0;
}

/*
 * pentatopes, p = 2
*/
template<>
void
Lagrange<Simplex,3,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  phi[0] =  s*-3.0-t*3.0-u*3.0+s*t*4.0+s*u*4.0+t*u*4.0+(s*s)*2.0+(t*t)*2.0+(u*u)*2.0+1.0;
  phi[1] =  -s+(s*s)*2.0;
  phi[2] =  -t+(t*t)*2.0;
  phi[3] =  -u+(u*u)*2.0;
  phi[4] =  t*u*4.0;
  phi[5] =  s*u*4.0;
  phi[6] =  s*t*4.0;
  phi[7] =  t*4.0-s*t*4.0-t*u*4.0-(t*t)*4.0;
  phi[8] =  u*4.0-s*u*4.0-t*u*4.0-(u*u)*4.0;
  phi[9] =  s*4.0-s*t*4.0-s*u*4.0-(s*s)*4.0;
}

template<>
void
Lagrange<Simplex,3,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t* phis = phi;
  phis[0] =  s*4.0+t*4.0+u*4.0-3.0;
  phis[1] =  s*4.0-1.0;
  phis[2] =  0.0;
  phis[3] =  0.0;
  phis[4] =  0.0;
  phis[5] =  u*4.0;
  phis[6] =  t*4.0;
  phis[7] =  t*-4.0;
  phis[8] =  u*-4.0;
  phis[9] =  s*-8.0-t*4.0-u*4.0+4.0;

  real_t* phit = phi + 10;
  phit[0] =  s*4.0+t*4.0+u*4.0-3.0;
  phit[1] =  0.0;
  phit[2] =  t*4.0-1.0;
  phit[3] =  0.0;
  phit[4] =  u*4.0;
  phit[5] =  0.0;
  phit[6] =  s*4.0;
  phit[7] =  s*-4.0-t*8.0-u*4.0+4.0;
  phit[8] =  u*-4.0;
  phit[9] =  s*-4.0;

  real_t* phiu = phi + 20;
  phiu[0] =  s*4.0+t*4.0+u*4.0-3.0;
  phiu[1] =  0.0;
  phiu[2] =  0.0;
  phiu[3] =  u*4.0-1.0;
  phiu[4] =  t*4.0;
  phiu[5] =  s*4.0;
  phiu[6] =  0.0;
  phiu[7] =  t*-4.0;
  phiu[8] =  s*-4.0-t*4.0-u*8.0+4.0;
  phiu[9] =  s*-4.0;

}

template<>
void
Lagrange<Simplex,3,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0] =  4.0;
  phiss[1] =  4.0;
  phiss[2] =  0.0;
  phiss[3] =  0.0;
  phiss[4] =  0.0;
  phiss[5] =  0.0;
  phiss[6] =  0.0;
  phiss[7] =  0.0;
  phiss[8] =  0.0;
  phiss[9] =  -8.0;

  // phist
  real_t *phist = phi + 10;
  phist[0] =  4.0;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  0.0;
  phist[4] =  0.0;
  phist[5] =  0.0;
  phist[6] =  4.0;
  phist[7] =  -4.0;
  phist[8] =  0.0;
  phist[9] =  -4.0;

  // phitt
  real_t *phitt = phi + 20;
  phitt[0] =  4.0;
  phitt[1] =  0.0;
  phitt[2] =  4.0;
  phitt[3] =  0.0;
  phitt[4] =  0.0;
  phitt[5] =  0.0;
  phitt[6] =  0.0;
  phitt[7] =  -8.0;
  phitt[8] =  0.0;
  phitt[9] =  0.0;

  // phisu
  real_t* phisu = phi + 30;
  phisu[0] =  4.0;
  phisu[1] =  0.0;
  phisu[2] =  0.0;
  phisu[3] =  0.0;
  phisu[4] =  0.0;
  phisu[5] =  4.0;
  phisu[6] =  0.0;
  phisu[7] =  0.0;
  phisu[8] =  -4.0;
  phisu[9] =  -4.0;

  // phitu
  real_t *phitu = phi + 40;
  phitu[0] =  4.0;
  phitu[1] =  0.0;
  phitu[2] =  0.0;
  phitu[3] =  0.0;
  phitu[4] =  4.0;
  phitu[5] =  0.0;
  phitu[6] =  0.0;
  phitu[7] =  -4.0;
  phitu[8] =  -4.0;
  phitu[9] =  0.0;

  // phiuu
  real_t *phiuu = phi + 50;
  phiuu[0] =  4.0;
  phiuu[1] =  0.0;
  phiuu[2] =  0.0;
  phiuu[3] =  4.0;
  phiuu[4] =  0.0;
  phiuu[5] =  0.0;
  phiuu[6] =  0.0;
  phiuu[7] =  0.0;
  phiuu[8] =  -8.0;
  phiuu[9] =  0.0;
}

/*
 * tetrahedron, p = 3
*/
template<>
void
Lagrange<Simplex,3,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)-u*(1.1E1/2.0)+s*t*1.8E1+s*u*1.8E1+t*u*1.8E1-s*(t*t)*(2.7E1/2.0)-
             (s*s)*t*(2.7E1/2.0)-s*(u*u)*(2.7E1/2.0)-(s*s)*u*(2.7E1/2.0)-t*(u*u)*(2.7E1/2.0)-(t*t)*u*(2.7E1/2.0)+(s*
             s)*9.0-(s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+(u*u)*9.0-(u*u*u)*(9.0/2.0)-s*t*u*2.7E1+1.0;
  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);
  phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);
  phi[3] =  u-(u*u)*(9.0/2.0)+(u*u*u)*(9.0/2.0);
  phi[4] =  t*u*(-9.0/2.0)+(t*t)*u*(2.7E1/2.0);
  phi[5] =  t*u*(-9.0/2.0)+t*(u*u)*(2.7E1/2.0);
  phi[6] =  s*u*(-9.0/2.0)+s*(u*u)*(2.7E1/2.0);
  phi[7] =  s*u*(-9.0/2.0)+(s*s)*u*(2.7E1/2.0);
  phi[8] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);
  phi[9] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);
  phi[10] =  t*(-9.0/2.0)+s*t*(9.0/2.0)+t*u*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)-(t*t)*u*(2.7E1/2.0)+(t*t)*1.8E1-
             (t*t*t)*(2.7E1/2.0);
  phi[11] =  t*9.0-s*t*(4.5E1/2.0)-t*u*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)+t*(u*u)*(2.7E1/2.0)+
             (t*t)*u*2.7E1-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0)+s*t*u*2.7E1;
  phi[12] =  u*9.0-s*u*(4.5E1/2.0)-t*u*(4.5E1/2.0)+s*(u*u)*2.7E1+(s*s)*u*(2.7E1/2.0)+t*(u*u)*2.7E1+(t*t)*
             u*(2.7E1/2.0)-(u*u)*(4.5E1/2.0)+(u*u*u)*(2.7E1/2.0)+s*t*u*2.7E1;
  phi[13] =  u*(-9.0/2.0)+s*u*(9.0/2.0)+t*u*(9.0/2.0)-s*(u*u)*(2.7E1/2.0)-t*(u*u)*(2.7E1/2.0)+(u*u)*1.8E1-
             (u*u*u)*(2.7E1/2.0);
  phi[14] =  s*9.0-s*t*(4.5E1/2.0)-s*u*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1+s*(u*u)*(2.7E1/2.0)+
             (s*s)*u*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0)+s*t*u*2.7E1;
  phi[15] =  s*(-9.0/2.0)+s*t*(9.0/2.0)+s*u*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)-(s*s)*u*(2.7E1/2.0)+(s*s)*1.8E1-
             (s*s*s)*(2.7E1/2.0);
  phi[16] =  s*t*u*2.7E1;
  phi[17] =  t*u*2.7E1-t*(u*u)*2.7E1-(t*t)*u*2.7E1-s*t*u*2.7E1;
  phi[18] =  s*u*2.7E1-s*(u*u)*2.7E1-(s*s)*u*2.7E1-s*t*u*2.7E1;
  phi[19] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1-s*t*u*2.7E1;
}

template<>
void
Lagrange<Simplex,3,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t* phis = phi;
  phis[0] =  s*1.8E1+t*1.8E1+u*1.8E1-s*t*2.7E1-s*u*2.7E1-t*u*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-
             (u*u)*(2.7E1/2.0)-1.1E1/2.0;
  phis[1] =  s*-9.0+(s*s)*(2.7E1/2.0)+1.0;
  phis[2] =  0.0;
  phis[3] =  0.0;
  phis[4] =  0.0;
  phis[5] =  0.0;
  phis[6] =  u*(-9.0/2.0)+(u*u)*(2.7E1/2.0);
  phis[7] =  u*(-9.0/2.0)+s*u*2.7E1;
  phis[8] =  t*(-9.0/2.0)+s*t*2.7E1;
  phis[9] =  t*(-9.0/2.0)+(t*t)*(2.7E1/2.0);
  phis[10] =  t*(9.0/2.0)-(t*t)*(2.7E1/2.0);
  phis[11] =  t*(-4.5E1/2.0)+s*t*2.7E1+t*u*2.7E1+(t*t)*2.7E1;
  phis[12] =  u*(-4.5E1/2.0)+s*u*2.7E1+t*u*2.7E1+(u*u)*2.7E1;
  phis[13] =  u*(9.0/2.0)-(u*u)*(2.7E1/2.0);
  phis[14] =  s*-4.5E1-t*(4.5E1/2.0)-u*(4.5E1/2.0)+s*t*5.4E1+s*u*5.4E1+t*u*2.7E1+(s*s)*(8.1E1/2.0)+(t*t)*
             (2.7E1/2.0)+(u*u)*(2.7E1/2.0)+9.0;
  phis[15] =  s*3.6E1+t*(9.0/2.0)+u*(9.0/2.0)-s*t*2.7E1-s*u*2.7E1-(s*s)*(8.1E1/2.0)-9.0/2.0;
  phis[16] =  t*u*2.7E1;
  phis[17] =  t*u*-2.7E1;
  phis[18] =  u*2.7E1-s*u*5.4E1-t*u*2.7E1-(u*u)*2.7E1;
  phis[19] =  t*2.7E1-s*t*5.4E1-t*u*2.7E1-(t*t)*2.7E1;

  real_t* phit = phi + 20;
  phit[0] =  s*1.8E1+t*1.8E1+u*1.8E1-s*t*2.7E1-s*u*2.7E1-t*u*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-
             (u*u)*(2.7E1/2.0)-1.1E1/2.0;
  phit[1] =  0.0;
  phit[2] =  t*-9.0+(t*t)*(2.7E1/2.0)+1.0;
  phit[3] =  0.0;
  phit[4] =  u*(-9.0/2.0)+t*u*2.7E1;
  phit[5] =  u*(-9.0/2.0)+(u*u)*(2.7E1/2.0);
  phit[6] =  0.0;
  phit[7] =  0.0;
  phit[8] =  s*(-9.0/2.0)+(s*s)*(2.7E1/2.0);
  phit[9] =  s*(-9.0/2.0)+s*t*2.7E1;
  phit[10] =  s*(9.0/2.0)+t*3.6E1+u*(9.0/2.0)-s*t*2.7E1-t*u*2.7E1-(t*t)*(8.1E1/2.0)-9.0/2.0;
  phit[11] =  s*(-4.5E1/2.0)-t*4.5E1-u*(4.5E1/2.0)+s*t*5.4E1+s*u*2.7E1+t*u*5.4E1+(s*s)*(2.7E1/2.0)+(t*t)*
             (8.1E1/2.0)+(u*u)*(2.7E1/2.0)+9.0;
  phit[12] =  u*(-4.5E1/2.0)+s*u*2.7E1+t*u*2.7E1+(u*u)*2.7E1;
  phit[13] =  u*(9.0/2.0)-(u*u)*(2.7E1/2.0);
  phit[14] =  s*(-4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+(s*s)*2.7E1;
  phit[15] =  s*(9.0/2.0)-(s*s)*(2.7E1/2.0);
  phit[16] =  s*u*2.7E1;
  phit[17] =  u*2.7E1-s*u*2.7E1-t*u*5.4E1-(u*u)*2.7E1;
  phit[18] =  s*u*-2.7E1;
  phit[19] =  s*2.7E1-s*t*5.4E1-s*u*2.7E1-(s*s)*2.7E1;

  real_t* phiu = phi + 40;
  phiu[0] =  s*1.8E1+t*1.8E1+u*1.8E1-s*t*2.7E1-s*u*2.7E1-t*u*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-
             (u*u)*(2.7E1/2.0)-1.1E1/2.0;
  phiu[1] =  0.0;
  phiu[2] =  0.0;
  phiu[3] =  u*-9.0+(u*u)*(2.7E1/2.0)+1.0;
  phiu[4] =  t*(-9.0/2.0)+(t*t)*(2.7E1/2.0);
  phiu[5] =  t*(-9.0/2.0)+t*u*2.7E1;
  phiu[6] =  s*(-9.0/2.0)+s*u*2.7E1;
  phiu[7] =  s*(-9.0/2.0)+(s*s)*(2.7E1/2.0);
  phiu[8] =  0.0;
  phiu[9] =  0.0;
  phiu[10] =  t*(9.0/2.0)-(t*t)*(2.7E1/2.0);
  phiu[11] =  t*(-4.5E1/2.0)+s*t*2.7E1+t*u*2.7E1+(t*t)*2.7E1;
  phiu[12] =  s*(-4.5E1/2.0)-t*(4.5E1/2.0)-u*4.5E1+s*t*2.7E1+s*u*5.4E1+t*u*5.4E1+(s*s)*(2.7E1/2.0)+(t*t)*
             (2.7E1/2.0)+(u*u)*(8.1E1/2.0)+9.0;
  phiu[13] =  s*(9.0/2.0)+t*(9.0/2.0)+u*3.6E1-s*u*2.7E1-t*u*2.7E1-(u*u)*(8.1E1/2.0)-9.0/2.0;
  phiu[14] =  s*(-4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+(s*s)*2.7E1;
  phiu[15] =  s*(9.0/2.0)-(s*s)*(2.7E1/2.0);
  phiu[16] =  s*t*2.7E1;
  phiu[17] =  t*2.7E1-s*t*2.7E1-t*u*5.4E1-(t*t)*2.7E1;
  phiu[18] =  s*2.7E1-s*t*2.7E1-s*u*5.4E1-(s*s)*2.7E1;
  phiu[19] =  s*t*-2.7E1;
}

template<>
void
Lagrange<Simplex,3,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  // phiss
  real_t *phiss = phi;
  phiss[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phiss[1] =  s*2.7E1-9.0;
  phiss[2] =  0.0;
  phiss[3] =  0.0;
  phiss[4] =  0.0;
  phiss[5] =  0.0;
  phiss[6] =  0.0;
  phiss[7] =  u*2.7E1;
  phiss[8] =  t*2.7E1;
  phiss[9] =  0.0;
  phiss[10] =  0.0;
  phiss[11] =  t*2.7E1;
  phiss[12] =  u*2.7E1;
  phiss[13] =  0.0;
  phiss[14] =  s*8.1E1+t*5.4E1+u*5.4E1-4.5E1;
  phiss[15] =  s*-8.1E1-t*2.7E1-u*2.7E1+3.6E1;
  phiss[16] =  0.0;
  phiss[17] =  0.0;
  phiss[18] =  u*-5.4E1;
  phiss[19] =  t*-5.4E1;

  // phist
  real_t *phist = phi + 20;
  phist[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  0.0;
  phist[4] =  0.0;
  phist[5] =  0.0;
  phist[6] =  0.0;
  phist[7] =  0.0;
  phist[8] =  s*2.7E1-9.0/2.0;
  phist[9] =  t*2.7E1-9.0/2.0;
  phist[10] =  t*-2.7E1+9.0/2.0;
  phist[11] =  s*2.7E1+t*5.4E1+u*2.7E1-4.5E1/2.0;
  phist[12] =  u*2.7E1;
  phist[13] =  0.0;
  phist[14] =  s*5.4E1+t*2.7E1+u*2.7E1-4.5E1/2.0;
  phist[15] =  s*-2.7E1+9.0/2.0;
  phist[16] =  u*2.7E1;
  phist[17] =  u*-2.7E1;
  phist[18] =  u*-2.7E1;
  phist[19] =  s*-5.4E1-t*5.4E1-u*2.7E1+2.7E1;

  // phitt
  real_t* phitt = phi + 40;
  phitt[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phitt[1] =  0.0;
  phitt[2] =  t*2.7E1-9.0;
  phitt[3] =  0.0;
  phitt[4] =  u*2.7E1;
  phitt[5] =  0.0;
  phitt[6] =  0.0;
  phitt[7] =  0.0;
  phitt[8] =  0.0;
  phitt[9] =  s*2.7E1;
  phitt[10] =  s*-2.7E1-t*8.1E1-u*2.7E1+3.6E1;
  phitt[11] =  s*5.4E1+t*8.1E1+u*5.4E1-4.5E1;
  phitt[12] =  u*2.7E1;
  phitt[13] =  0.0;
  phitt[14] =  s*2.7E1;
  phitt[15] =  0.0;
  phitt[16] =  0.0;
  phitt[17] =  u*-5.4E1;
  phitt[18] =  0.0;
  phitt[19] =  s*-5.4E1;

  // phisu
  real_t* phisu = phi + 60;
  phisu[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phisu[1] =  0.0;
  phisu[2] =  0.0;
  phisu[3] =  0.0;
  phisu[4] =  0.0;
  phisu[5] =  0.0;
  phisu[6] =  u*2.7E1-9.0/2.0;
  phisu[7] =  s*2.7E1-9.0/2.0;
  phisu[8] =  0.0;
  phisu[9] =  0.0;
  phisu[10] =  0.0;
  phisu[11] =  t*2.7E1;
  phisu[12] =  s*2.7E1+t*2.7E1+u*5.4E1-4.5E1/2.0;
  phisu[13] =  u*-2.7E1+9.0/2.0;
  phisu[14] =  s*5.4E1+t*2.7E1+u*2.7E1-4.5E1/2.0;
  phisu[15] =  s*-2.7E1+9.0/2.0;
  phisu[16] =  t*2.7E1;
  phisu[17] =  t*-2.7E1;
  phisu[18] =  s*-5.4E1-t*2.7E1-u*5.4E1+2.7E1;
  phisu[19] =  t*-2.7E1;

  // phitu
  real_t *phitu = phi + 80;
  phitu[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phitu[1] =  0.0;
  phitu[2] =  0.0;
  phitu[3] =  0.0;
  phitu[4] =  t*2.7E1-9.0/2.0;
  phitu[5] =  u*2.7E1-9.0/2.0;
  phitu[6] =  0.0;
  phitu[7] =  0.0;
  phitu[8] =  0.0;
  phitu[9] =  0.0;
  phitu[10] =  t*-2.7E1+9.0/2.0;
  phitu[11] =  s*2.7E1+t*5.4E1+u*2.7E1-4.5E1/2.0;
  phitu[12] =  s*2.7E1+t*2.7E1+u*5.4E1-4.5E1/2.0;
  phitu[13] =  u*-2.7E1+9.0/2.0;
  phitu[14] =  s*2.7E1;
  phitu[15] =  0.0;
  phitu[16] =  s*2.7E1;
  phitu[17] =  s*-2.7E1-t*5.4E1-u*5.4E1+2.7E1;
  phitu[18] =  s*-2.7E1;
  phitu[19] =  s*-2.7E1;

  // phiuu
  real_t *phiuu = phi + 100;
  phiuu[0] =  s*-2.7E1-t*2.7E1-u*2.7E1+1.8E1;
  phiuu[1] =  0.0;
  phiuu[2] =  0.0;
  phiuu[3] =  u*2.7E1-9.0;
  phiuu[4] =  0.0;
  phiuu[5] =  t*2.7E1;
  phiuu[6] =  s*2.7E1;
  phiuu[7] =  0.0;
  phiuu[8] =  0.0;
  phiuu[9] =  0.0;
  phiuu[10] =  0.0;
  phiuu[11] =  t*2.7E1;
  phiuu[12] =  s*5.4E1+t*5.4E1+u*8.1E1-4.5E1;
  phiuu[13] =  s*-2.7E1-t*2.7E1-u*8.1E1+3.6E1;
  phiuu[14] =  s*2.7E1;
  phiuu[15] =  0.0;
  phiuu[16] =  0.0;
  phiuu[17] =  t*-5.4E1;
  phiuu[18] =  s*-5.4E1;
  phiuu[19] =  0.0;
}

/*
 * tetrahedra, p = 4
*/
template<>
void
Lagrange<Simplex,3,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  // phi
  phi[0] =  (s*s)*(t*t)*6.4E1+(s*s)*(u*u)*6.4E1+(t*t)*(u*u)*6.4E1-ONETHIRD*s*2.5E1-ONETHIRD*t*2.5E1-ONETHIRD*
             u*2.5E1+ONETHIRD*(s*s)*7.0E1-ONETHIRD*(s*s*s)*8.0E1+ONETHIRD*(s*s*s*s)*3.2E1+ONETHIRD*(t*t)*7.0E1-ONETHIRD*
             (t*t*t)*8.0E1+ONETHIRD*(t*t*t*t)*3.2E1+ONETHIRD*(u*u)*7.0E1-ONETHIRD*(u*u*u)*8.0E1+ONETHIRD*(u*u*u*u)*
             3.2E1-s*(t*t)*8.0E1-(s*s)*t*8.0E1-s*(u*u)*8.0E1-(s*s)*u*8.0E1-t*(u*u)*8.0E1-(t*t)*u*8.0E1+s*t*(u*u)*1.28E2+
             s*(t*t)*u*1.28E2+(s*s)*t*u*1.28E2+ONETHIRD*s*t*1.4E2+ONETHIRD*s*u*1.4E2+ONETHIRD*t*u*1.4E2-s*t*u*1.6E2+
             ONETHIRD*s*(t*t*t)*1.28E2+ONETHIRD*(s*s*s)*t*1.28E2+ONETHIRD*s*(u*u*u)*1.28E2+ONETHIRD*(s*s*s)*u*1.28E2+
             ONETHIRD*t*(u*u*u)*1.28E2+ONETHIRD*(t*t*t)*u*1.28E2+1.0;
  phi[1] =  -s+ONETHIRD*(s*s)*2.2E1+ONETHIRD*(s*s*s*s)*3.2E1-(s*s*s)*1.6E1;
  phi[2] =  -t+ONETHIRD*(t*t)*2.2E1+ONETHIRD*(t*t*t*t)*3.2E1-(t*t*t)*1.6E1;
  phi[3] =  -u+ONETHIRD*(u*u)*2.2E1+ONETHIRD*(u*u*u*u)*3.2E1-(u*u*u)*1.6E1;
  phi[4] =  (t*t)*u*-3.2E1+ONETHIRD*t*u*1.6E1+ONETHIRD*(t*t*t)*u*1.28E2;
  phi[5] =  (t*t)*(u*u)*6.4E1+t*u*4.0-t*(u*u)*1.6E1-(t*t)*u*1.6E1;
  phi[6] =  t*(u*u)*-3.2E1+ONETHIRD*t*u*1.6E1+ONETHIRD*t*(u*u*u)*1.28E2;
  phi[7] =  s*(u*u)*-3.2E1+ONETHIRD*s*u*1.6E1+ONETHIRD*s*(u*u*u)*1.28E2;
  phi[8] =  (s*s)*(u*u)*6.4E1+s*u*4.0-s*(u*u)*1.6E1-(s*s)*u*1.6E1;
  phi[9] =  (s*s)*u*-3.2E1+ONETHIRD*s*u*1.6E1+ONETHIRD*(s*s*s)*u*1.28E2;
  phi[10] =  (s*s)*t*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*(s*s*s)*t*1.28E2;
  phi[11] =  (s*s)*(t*t)*6.4E1+s*t*4.0-s*(t*t)*1.6E1-(s*s)*t*1.6E1;
  phi[12] =  s*(t*t)*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*s*(t*t*t)*1.28E2;
  phi[13] =  ONETHIRD*t*1.6E1-ONETHIRD*(t*t)*1.12E2+ONETHIRD*(t*t*t)*2.24E2-ONETHIRD*(t*t*t*t)*1.28E2+s*
             (t*t)*3.2E1+(t*t)*u*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*t*u*1.6E1-ONETHIRD*s*(t*t*t)*1.28E2-ONETHIRD*(t*
             t*t)*u*1.28E2;
  phi[14] =  t*-1.2E1+(s*s)*(t*t)*6.4E1+(t*t)*(u*u)*6.4E1+s*t*2.8E1+t*u*2.8E1-s*(t*t)*1.44E2-(s*s)*t*1.6E1+
             s*(t*t*t)*1.28E2-t*(u*u)*1.6E1-(t*t)*u*1.44E2+(t*t*t)*u*1.28E2+(t*t)*7.6E1-(t*t*t)*1.28E2+(t*t*t*t)*6.4E1+
             s*(t*t)*u*1.28E2-s*t*u*3.2E1;
  phi[15] =  t*1.6E1-(s*s)*(t*t)*1.28E2-(t*t)*(u*u)*1.28E2-ONETHIRD*(t*t)*2.08E2-ONETHIRD*(t*t*t*t)*1.28E2+
             s*(t*t)*1.92E2+(s*s)*t*9.6E1-s*(t*t*t)*1.28E2+t*(u*u)*9.6E1+(t*t)*u*1.92E2-(t*t*t)*u*1.28E2+(t*t*t)*9.6E1-
             s*t*(u*u)*1.28E2-s*(t*t)*u*2.56E2-(s*s)*t*u*1.28E2-ONETHIRD*s*t*2.08E2-ONETHIRD*t*u*2.08E2+s*t*u*1.92E2-
             ONETHIRD*(s*s*s)*t*1.28E2-ONETHIRD*t*(u*u*u)*1.28E2;
  phi[16] =  u*1.6E1-(s*s)*(u*u)*1.28E2-(t*t)*(u*u)*1.28E2-ONETHIRD*(u*u)*2.08E2-ONETHIRD*(u*u*u*u)*1.28E2+
             s*(u*u)*1.92E2+(s*s)*u*9.6E1-s*(u*u*u)*1.28E2+t*(u*u)*1.92E2+(t*t)*u*9.6E1-t*(u*u*u)*1.28E2+(u*u*u)*9.6E1-
             s*t*(u*u)*2.56E2-s*(t*t)*u*1.28E2-(s*s)*t*u*1.28E2-ONETHIRD*s*u*2.08E2-ONETHIRD*t*u*2.08E2+s*t*u*1.92E2-
             ONETHIRD*(s*s*s)*u*1.28E2-ONETHIRD*(t*t*t)*u*1.28E2;
  phi[17] =  u*-1.2E1+(s*s)*(u*u)*6.4E1+(t*t)*(u*u)*6.4E1+s*u*2.8E1+t*u*2.8E1-s*(u*u)*1.44E2-(s*s)*u*1.6E1+
             s*(u*u*u)*1.28E2-t*(u*u)*1.44E2-(t*t)*u*1.6E1+t*(u*u*u)*1.28E2+(u*u)*7.6E1-(u*u*u)*1.28E2+(u*u*u*u)*6.4E1+
             s*t*(u*u)*1.28E2-s*t*u*3.2E1;
  phi[18] =  ONETHIRD*u*1.6E1-ONETHIRD*(u*u)*1.12E2+ONETHIRD*(u*u*u)*2.24E2-ONETHIRD*(u*u*u*u)*1.28E2+s*
             (u*u)*3.2E1+t*(u*u)*3.2E1-ONETHIRD*s*u*1.6E1-ONETHIRD*t*u*1.6E1-ONETHIRD*s*(u*u*u)*1.28E2-ONETHIRD*t*
             (u*u*u)*1.28E2;
  phi[19] =  s*1.6E1-(s*s)*(t*t)*1.28E2-(s*s)*(u*u)*1.28E2-ONETHIRD*(s*s)*2.08E2-ONETHIRD*(s*s*s*s)*1.28E2+
             s*(t*t)*9.6E1+(s*s)*t*1.92E2-(s*s*s)*t*1.28E2+s*(u*u)*9.6E1+(s*s)*u*1.92E2-(s*s*s)*u*1.28E2+(s*s*s)*9.6E1-
             s*t*(u*u)*1.28E2-s*(t*t)*u*1.28E2-(s*s)*t*u*2.56E2-ONETHIRD*s*t*2.08E2-ONETHIRD*s*u*2.08E2+s*t*u*1.92E2-
             ONETHIRD*s*(t*t*t)*1.28E2-ONETHIRD*s*(u*u*u)*1.28E2;
  phi[20] =  s*-1.2E1+(s*s)*(t*t)*6.4E1+(s*s)*(u*u)*6.4E1+s*t*2.8E1+s*u*2.8E1-s*(t*t)*1.6E1-(s*s)*t*1.44E2+
             (s*s*s)*t*1.28E2-s*(u*u)*1.6E1-(s*s)*u*1.44E2+(s*s*s)*u*1.28E2+(s*s)*7.6E1-(s*s*s)*1.28E2+(s*s*s*s)*6.4E1+
             (s*s)*t*u*1.28E2-s*t*u*3.2E1;
  phi[21] =  ONETHIRD*s*1.6E1-ONETHIRD*(s*s)*1.12E2+ONETHIRD*(s*s*s)*2.24E2-ONETHIRD*(s*s*s*s)*1.28E2+(s*
             s)*t*3.2E1+(s*s)*u*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*s*u*1.6E1-ONETHIRD*(s*s*s)*t*1.28E2-ONETHIRD*(s*
             s*s)*u*1.28E2;
  phi[22] =  (s*s)*t*u*1.28E2-s*t*u*3.2E1;
  phi[23] =  s*(t*t)*u*1.28E2-s*t*u*3.2E1;
  phi[24] =  s*t*(u*u)*1.28E2-s*t*u*3.2E1;
  phi[25] =  (t*t)*(u*u)*2.56E2+t*u*9.6E1-t*(u*u)*2.24E2-(t*t)*u*2.24E2+t*(u*u*u)*1.28E2+(t*t*t)*u*1.28E2+
             s*t*(u*u)*2.56E2+s*(t*t)*u*2.56E2+(s*s)*t*u*1.28E2-s*t*u*2.24E2;
  phi[26] =  (t*t)*(u*u)*-1.28E2-t*u*3.2E1+t*(u*u)*1.6E2+(t*t)*u*3.2E1-t*(u*u*u)*1.28E2-s*t*(u*u)*1.28E2+
             s*t*u*3.2E1;
  phi[27] =  (t*t)*(u*u)*-1.28E2-t*u*3.2E1+t*(u*u)*3.2E1+(t*t)*u*1.6E2-(t*t*t)*u*1.28E2-s*(t*t)*u*1.28E2+
             s*t*u*3.2E1;
  phi[28] =  (s*s)*(u*u)*2.56E2+s*u*9.6E1-s*(u*u)*2.24E2-(s*s)*u*2.24E2+s*(u*u*u)*1.28E2+(s*s*s)*u*1.28E2+
             s*t*(u*u)*2.56E2+s*(t*t)*u*1.28E2+(s*s)*t*u*2.56E2-s*t*u*2.24E2;
  phi[29] =  (s*s)*(u*u)*-1.28E2-s*u*3.2E1+s*(u*u)*3.2E1+(s*s)*u*1.6E2-(s*s*s)*u*1.28E2-(s*s)*t*u*1.28E2+
             s*t*u*3.2E1;
  phi[30] =  (s*s)*(u*u)*-1.28E2-s*u*3.2E1+s*(u*u)*1.6E2+(s*s)*u*3.2E1-s*(u*u*u)*1.28E2-s*t*(u*u)*1.28E2+
             s*t*u*3.2E1;
  phi[31] =  (s*s)*(t*t)*2.56E2+s*t*9.6E1-s*(t*t)*2.24E2-(s*s)*t*2.24E2+s*(t*t*t)*1.28E2+(s*s*s)*t*1.28E2+
             s*t*(u*u)*1.28E2+s*(t*t)*u*2.56E2+(s*s)*t*u*2.56E2-s*t*u*2.24E2;
  phi[32] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*1.6E2+(s*s)*t*3.2E1-s*(t*t*t)*1.28E2-s*(t*t)*u*1.28E2+
             s*t*u*3.2E1;
  phi[33] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*3.2E1+(s*s)*t*1.6E2-(s*s*s)*t*1.28E2-(s*s)*t*u*1.28E2+
             s*t*u*3.2E1;
  phi[34] =  s*t*(u*u)*-2.56E2-s*(t*t)*u*2.56E2-(s*s)*t*u*2.56E2+s*t*u*2.56E2;
}

template<>
void
Lagrange<Simplex,3,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  // phis
  real_t* phis = phi;
  phis[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2+ONETHIRD*u*1.4E2-s*t*1.6E2-s*u*1.6E2-t*u*
           1.6E2-ONETHIRD*(s*s)*2.4E2+ONETHIRD*(s*s*s)*1.28E2+ONETHIRD*(t*t*t)*1.28E2+ONETHIRD*(u*u*u)*1.28E2+s*
           (t*t)*1.28E2+s*(u*u)*1.28E2+t*(u*u)*1.28E2+(t*t)*u*1.28E2-(t*t)*8.0E1-(u*u)*8.0E1+s*t*u*2.56E2+ONETHIRD*
           (s*s)*t*3.84E2+ONETHIRD*(s*s)*u*3.84E2;
  phis[1] =  ONETHIRD*s*4.4E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*4.8E1-1.0;
  phis[2] =  0.0;
  phis[3] =  0.0;
  phis[4] =  0.0;
  phis[5] =  0.0;
  phis[6] =  0.0;
  phis[7] =  ONETHIRD*u*1.6E1+ONETHIRD*(u*u*u)*1.28E2-(u*u)*3.2E1;
  phis[8] =  u*4.0-s*u*3.2E1+s*(u*u)*1.28E2-(u*u)*1.6E1;
  phis[9] =  ONETHIRD*u*1.6E1-s*u*6.4E1+ONETHIRD*(s*s)*u*3.84E2;
  phis[10] =  ONETHIRD*t*1.6E1-s*t*6.4E1+ONETHIRD*(s*s)*t*3.84E2;
  phis[11] =  t*4.0-s*t*3.2E1+s*(t*t)*1.28E2-(t*t)*1.6E1;
  phis[12] =  ONETHIRD*t*1.6E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*3.2E1;
  phis[13] =  ONETHIRD*t*-1.6E1-ONETHIRD*(t*t*t)*1.28E2+(t*t)*3.2E1;
  phis[14] =  t*2.8E1-s*t*3.2E1-t*u*3.2E1+s*(t*t)*1.28E2+(t*t)*u*1.28E2-(t*t)*1.44E2+(t*t*t)*1.28E2;
  phis[15] =  ONETHIRD*t*-2.08E2+s*t*1.92E2+t*u*1.92E2-s*(t*t)*2.56E2-t*(u*u)*1.28E2-(t*t)*u*2.56E2+(t*
           t)*1.92E2-(t*t*t)*1.28E2-s*t*u*2.56E2-ONETHIRD*(s*s)*t*3.84E2;
  phis[16] =  ONETHIRD*u*-2.08E2+s*u*1.92E2+t*u*1.92E2-s*(u*u)*2.56E2-t*(u*u)*2.56E2-(t*t)*u*1.28E2+(u*
           u)*1.92E2-(u*u*u)*1.28E2-s*t*u*2.56E2-ONETHIRD*(s*s)*u*3.84E2;
  phis[17] =  u*2.8E1-s*u*3.2E1-t*u*3.2E1+s*(u*u)*1.28E2+t*(u*u)*1.28E2-(u*u)*1.44E2+(u*u*u)*1.28E2;
  phis[18] =  ONETHIRD*u*-1.6E1-ONETHIRD*(u*u*u)*1.28E2+(u*u)*3.2E1;
  phis[19] =  ONETHIRD*s*-4.16E2-ONETHIRD*t*2.08E2-ONETHIRD*u*2.08E2+s*t*3.84E2+s*u*3.84E2+t*u*1.92E2-ONETHIRD*
           (s*s*s)*5.12E2-ONETHIRD*(t*t*t)*1.28E2-ONETHIRD*(u*u*u)*1.28E2-s*(t*t)*2.56E2-(s*s)*t*3.84E2-s*(u*u)*
           2.56E2-(s*s)*u*3.84E2-t*(u*u)*1.28E2-(t*t)*u*1.28E2+(s*s)*2.88E2+(t*t)*9.6E1+(u*u)*9.6E1-s*t*u*5.12E2+
           1.6E1;
  phis[20] =  s*1.52E2+t*2.8E1+u*2.8E1-s*t*2.88E2-s*u*2.88E2-t*u*3.2E1+s*(t*t)*1.28E2+(s*s)*t*3.84E2+s*
           (u*u)*1.28E2+(s*s)*u*3.84E2-(s*s)*3.84E2+(s*s*s)*2.56E2-(t*t)*1.6E1-(u*u)*1.6E1+s*t*u*2.56E2-1.2E1;
  phis[21] =  ONETHIRD*1.6E1-ONETHIRD*s*2.24E2-ONETHIRD*t*1.6E1-ONETHIRD*u*1.6E1+s*t*6.4E1+s*u*6.4E1+ONETHIRD*
           (s*s)*6.72E2-ONETHIRD*(s*s*s)*5.12E2-ONETHIRD*(s*s)*t*3.84E2-ONETHIRD*(s*s)*u*3.84E2;
  phis[22] =  t*u*-3.2E1+s*t*u*2.56E2;
  phis[23] =  t*u*-3.2E1+(t*t)*u*1.28E2;
  phis[24] =  t*u*-3.2E1+t*(u*u)*1.28E2;
  phis[25] =  t*u*-2.24E2+t*(u*u)*2.56E2+(t*t)*u*2.56E2+s*t*u*2.56E2;
  phis[26] =  t*u*3.2E1-t*(u*u)*1.28E2;
  phis[27] =  t*u*3.2E1-(t*t)*u*1.28E2;
  phis[28] =  u*9.6E1-s*u*4.48E2-t*u*2.24E2+s*(u*u)*5.12E2+(s*s)*u*3.84E2+t*(u*u)*2.56E2+(t*t)*u*1.28E2-
           (u*u)*2.24E2+(u*u*u)*1.28E2+s*t*u*5.12E2;
  phis[29] =  u*-3.2E1+s*u*3.2E2+t*u*3.2E1-s*(u*u)*2.56E2-(s*s)*u*3.84E2+(u*u)*3.2E1-s*t*u*2.56E2;
  phis[30] =  u*-3.2E1+s*u*6.4E1+t*u*3.2E1-s*(u*u)*2.56E2-t*(u*u)*1.28E2+(u*u)*1.6E2-(u*u*u)*1.28E2;
  phis[31] =  t*9.6E1-s*t*4.48E2-t*u*2.24E2+s*(t*t)*5.12E2+(s*s)*t*3.84E2+t*(u*u)*1.28E2+(t*t)*u*2.56E2-
           (t*t)*2.24E2+(t*t*t)*1.28E2+s*t*u*5.12E2;
  phis[32] =  t*-3.2E1+s*t*6.4E1+t*u*3.2E1-s*(t*t)*2.56E2-(t*t)*u*1.28E2+(t*t)*1.6E2-(t*t*t)*1.28E2;
  phis[33] =  t*-3.2E1+s*t*3.2E2+t*u*3.2E1-s*(t*t)*2.56E2-(s*s)*t*3.84E2+(t*t)*3.2E1-s*t*u*2.56E2;
  phis[34] =  t*u*2.56E2-t*(u*u)*2.56E2-(t*t)*u*2.56E2-s*t*u*5.12E2;

  // phit
  real_t *phit = phi + 35;
  phit[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2+ONETHIRD*u*1.4E2-s*t*1.6E2-s*u*1.6E2-t*u*
           1.6E2+ONETHIRD*(s*s*s)*1.28E2-ONETHIRD*(t*t)*2.4E2+ONETHIRD*(t*t*t)*1.28E2+ONETHIRD*(u*u*u)*1.28E2+(s*
           s)*t*1.28E2+s*(u*u)*1.28E2+(s*s)*u*1.28E2+t*(u*u)*1.28E2-(s*s)*8.0E1-(u*u)*8.0E1+s*t*u*2.56E2+ONETHIRD*
           s*(t*t)*3.84E2+ONETHIRD*(t*t)*u*3.84E2;
  phit[1] =  0.0;
  phit[2] =  ONETHIRD*t*4.4E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*4.8E1-1.0;
  phit[3] =  0.0;
  phit[4] =  ONETHIRD*u*1.6E1-t*u*6.4E1+ONETHIRD*(t*t)*u*3.84E2;
  phit[5] =  u*4.0-t*u*3.2E1+t*(u*u)*1.28E2-(u*u)*1.6E1;
  phit[6] =  ONETHIRD*u*1.6E1+ONETHIRD*(u*u*u)*1.28E2-(u*u)*3.2E1;
  phit[7] =  0.0;
  phit[8] =  0.0;
  phit[9] =  0.0;
  phit[10] =  ONETHIRD*s*1.6E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*3.2E1;
  phit[11] =  s*4.0-s*t*3.2E1+(s*s)*t*1.28E2-(s*s)*1.6E1;
  phit[12] =  ONETHIRD*s*1.6E1-s*t*6.4E1+ONETHIRD*s*(t*t)*3.84E2;
  phit[13] =  ONETHIRD*1.6E1-ONETHIRD*s*1.6E1-ONETHIRD*t*2.24E2-ONETHIRD*u*1.6E1+s*t*6.4E1+t*u*6.4E1+ONETHIRD*
           (t*t)*6.72E2-ONETHIRD*(t*t*t)*5.12E2-ONETHIRD*s*(t*t)*3.84E2-ONETHIRD*(t*t)*u*3.84E2;
  phit[14] =  s*2.8E1+t*1.52E2+u*2.8E1-s*t*2.88E2-s*u*3.2E1-t*u*2.88E2+s*(t*t)*3.84E2+(s*s)*t*1.28E2+t*
           (u*u)*1.28E2+(t*t)*u*3.84E2-(s*s)*1.6E1-(t*t)*3.84E2+(t*t*t)*2.56E2-(u*u)*1.6E1+s*t*u*2.56E2-1.2E1;
  phit[15] =  ONETHIRD*s*-2.08E2-ONETHIRD*t*4.16E2-ONETHIRD*u*2.08E2+s*t*3.84E2+s*u*1.92E2+t*u*3.84E2-ONETHIRD*
           (s*s*s)*1.28E2-ONETHIRD*(t*t*t)*5.12E2-ONETHIRD*(u*u*u)*1.28E2-s*(t*t)*3.84E2-(s*s)*t*2.56E2-s*(u*u)*
           1.28E2-(s*s)*u*1.28E2-t*(u*u)*2.56E2-(t*t)*u*3.84E2+(s*s)*9.6E1+(t*t)*2.88E2+(u*u)*9.6E1-s*t*u*5.12E2+
           1.6E1;
  phit[16] =  ONETHIRD*u*-2.08E2+s*u*1.92E2+t*u*1.92E2-s*(u*u)*2.56E2-(s*s)*u*1.28E2-t*(u*u)*2.56E2+(u*
           u)*1.92E2-(u*u*u)*1.28E2-s*t*u*2.56E2-ONETHIRD*(t*t)*u*3.84E2;
  phit[17] =  u*2.8E1-s*u*3.2E1-t*u*3.2E1+s*(u*u)*1.28E2+t*(u*u)*1.28E2-(u*u)*1.44E2+(u*u*u)*1.28E2;
  phit[18] =  ONETHIRD*u*-1.6E1-ONETHIRD*(u*u*u)*1.28E2+(u*u)*3.2E1;
  phit[19] =  ONETHIRD*s*-2.08E2+s*t*1.92E2+s*u*1.92E2-(s*s)*t*2.56E2-s*(u*u)*1.28E2-(s*s)*u*2.56E2+(s*
           s)*1.92E2-(s*s*s)*1.28E2-s*t*u*2.56E2-ONETHIRD*s*(t*t)*3.84E2;
  phit[20] =  s*2.8E1-s*t*3.2E1-s*u*3.2E1+(s*s)*t*1.28E2+(s*s)*u*1.28E2-(s*s)*1.44E2+(s*s*s)*1.28E2;
  phit[21] =  ONETHIRD*s*-1.6E1-ONETHIRD*(s*s*s)*1.28E2+(s*s)*3.2E1;
  phit[22] =  s*u*-3.2E1+(s*s)*u*1.28E2;
  phit[23] =  s*u*-3.2E1+s*t*u*2.56E2;
  phit[24] =  s*u*-3.2E1+s*(u*u)*1.28E2;
  phit[25] =  u*9.6E1-s*u*2.24E2-t*u*4.48E2+s*(u*u)*2.56E2+(s*s)*u*1.28E2+t*(u*u)*5.12E2+(t*t)*u*3.84E2-
           (u*u)*2.24E2+(u*u*u)*1.28E2+s*t*u*5.12E2;
  phit[26] =  u*-3.2E1+s*u*3.2E1+t*u*6.4E1-s*(u*u)*1.28E2-t*(u*u)*2.56E2+(u*u)*1.6E2-(u*u*u)*1.28E2;
  phit[27] =  u*-3.2E1+s*u*3.2E1+t*u*3.2E2-t*(u*u)*2.56E2-(t*t)*u*3.84E2+(u*u)*3.2E1-s*t*u*2.56E2;
  phit[28] =  s*u*-2.24E2+s*(u*u)*2.56E2+(s*s)*u*2.56E2+s*t*u*2.56E2;
  phit[29] =  s*u*3.2E1-(s*s)*u*1.28E2;
  phit[30] =  s*u*3.2E1-s*(u*u)*1.28E2;
  phit[31] =  s*9.6E1-s*t*4.48E2-s*u*2.24E2+s*(t*t)*3.84E2+(s*s)*t*5.12E2+s*(u*u)*1.28E2+(s*s)*u*2.56E2-
           (s*s)*2.24E2+(s*s*s)*1.28E2+s*t*u*5.12E2;
  phit[32] =  s*-3.2E1+s*t*3.2E2+s*u*3.2E1-s*(t*t)*3.84E2-(s*s)*t*2.56E2+(s*s)*3.2E1-s*t*u*2.56E2;
  phit[33] =  s*-3.2E1+s*t*6.4E1+s*u*3.2E1-(s*s)*t*2.56E2-(s*s)*u*1.28E2+(s*s)*1.6E2-(s*s*s)*1.28E2;
  phit[34] =  s*u*2.56E2-s*(u*u)*2.56E2-(s*s)*u*2.56E2-s*t*u*5.12E2;

  // phiu
  real_t* phiu = phi + 70;
  phiu[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2+ONETHIRD*u*1.4E2-s*t*1.6E2-s*u*1.6E2-t*u*
           1.6E2+ONETHIRD*(s*s*s)*1.28E2+ONETHIRD*(t*t*t)*1.28E2-ONETHIRD*(u*u)*2.4E2+ONETHIRD*(u*u*u)*1.28E2+s*
           (t*t)*1.28E2+(s*s)*t*1.28E2+(s*s)*u*1.28E2+(t*t)*u*1.28E2-(s*s)*8.0E1-(t*t)*8.0E1+s*t*u*2.56E2+ONETHIRD*
           s*(u*u)*3.84E2+ONETHIRD*t*(u*u)*3.84E2;
  phiu[1] =  0.0;
  phiu[2] =  0.0;
  phiu[3] =  ONETHIRD*u*4.4E1+ONETHIRD*(u*u*u)*1.28E2-(u*u)*4.8E1-1.0;
  phiu[4] =  ONETHIRD*t*1.6E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*3.2E1;
  phiu[5] =  t*4.0-t*u*3.2E1+(t*t)*u*1.28E2-(t*t)*1.6E1;
  phiu[6] =  ONETHIRD*t*1.6E1-t*u*6.4E1+ONETHIRD*t*(u*u)*3.84E2;
  phiu[7] =  ONETHIRD*s*1.6E1-s*u*6.4E1+ONETHIRD*s*(u*u)*3.84E2;
  phiu[8] =  s*4.0-s*u*3.2E1+(s*s)*u*1.28E2-(s*s)*1.6E1;
  phiu[9] =  ONETHIRD*s*1.6E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*3.2E1;
  phiu[10] =  0.0;
  phiu[11] =  0.0;
  phiu[12] =  0.0;
  phiu[13] =  ONETHIRD*t*-1.6E1-ONETHIRD*(t*t*t)*1.28E2+(t*t)*3.2E1;
  phiu[14] =  t*2.8E1-s*t*3.2E1-t*u*3.2E1+s*(t*t)*1.28E2+(t*t)*u*1.28E2-(t*t)*1.44E2+(t*t*t)*1.28E2;
  phiu[15] =  ONETHIRD*t*-2.08E2+s*t*1.92E2+t*u*1.92E2-s*(t*t)*2.56E2-(s*s)*t*1.28E2-(t*t)*u*2.56E2+(t*
           t)*1.92E2-(t*t*t)*1.28E2-s*t*u*2.56E2-ONETHIRD*t*(u*u)*3.84E2;
  phiu[16] =  ONETHIRD*s*-2.08E2-ONETHIRD*t*2.08E2-ONETHIRD*u*4.16E2+s*t*1.92E2+s*u*3.84E2+t*u*3.84E2-ONETHIRD*
           (s*s*s)*1.28E2-ONETHIRD*(t*t*t)*1.28E2-ONETHIRD*(u*u*u)*5.12E2-s*(t*t)*1.28E2-(s*s)*t*1.28E2-s*(u*u)*
           3.84E2-(s*s)*u*2.56E2-t*(u*u)*3.84E2-(t*t)*u*2.56E2+(s*s)*9.6E1+(t*t)*9.6E1+(u*u)*2.88E2-s*t*u*5.12E2+
           1.6E1;
  phiu[17] =  s*2.8E1+t*2.8E1+u*1.52E2-s*t*3.2E1-s*u*2.88E2-t*u*2.88E2+s*(u*u)*3.84E2+(s*s)*u*1.28E2+t*
           (u*u)*3.84E2+(t*t)*u*1.28E2-(s*s)*1.6E1-(t*t)*1.6E1-(u*u)*3.84E2+(u*u*u)*2.56E2+s*t*u*2.56E2-1.2E1;
  phiu[18] =  ONETHIRD*1.6E1-ONETHIRD*s*1.6E1-ONETHIRD*t*1.6E1-ONETHIRD*u*2.24E2+s*u*6.4E1+t*u*6.4E1+ONETHIRD*
           (u*u)*6.72E2-ONETHIRD*(u*u*u)*5.12E2-ONETHIRD*s*(u*u)*3.84E2-ONETHIRD*t*(u*u)*3.84E2;
  phiu[19] =  ONETHIRD*s*-2.08E2+s*t*1.92E2+s*u*1.92E2-s*(t*t)*1.28E2-(s*s)*t*2.56E2-(s*s)*u*2.56E2+(s*
           s)*1.92E2-(s*s*s)*1.28E2-s*t*u*2.56E2-ONETHIRD*s*(u*u)*3.84E2;
  phiu[20] =  s*2.8E1-s*t*3.2E1-s*u*3.2E1+(s*s)*t*1.28E2+(s*s)*u*1.28E2-(s*s)*1.44E2+(s*s*s)*1.28E2;
  phiu[21] =  ONETHIRD*s*-1.6E1-ONETHIRD*(s*s*s)*1.28E2+(s*s)*3.2E1;
  phiu[22] =  s*t*-3.2E1+(s*s)*t*1.28E2;
  phiu[23] =  s*t*-3.2E1+s*(t*t)*1.28E2;
  phiu[24] =  s*t*-3.2E1+s*t*u*2.56E2;
  phiu[25] =  t*9.6E1-s*t*2.24E2-t*u*4.48E2+s*(t*t)*2.56E2+(s*s)*t*1.28E2+t*(u*u)*3.84E2+(t*t)*u*5.12E2-
           (t*t)*2.24E2+(t*t*t)*1.28E2+s*t*u*5.12E2;
  phiu[26] =  t*-3.2E1+s*t*3.2E1+t*u*3.2E2-t*(u*u)*3.84E2-(t*t)*u*2.56E2+(t*t)*3.2E1-s*t*u*2.56E2;
  phiu[27] =  t*-3.2E1+s*t*3.2E1+t*u*6.4E1-s*(t*t)*1.28E2-(t*t)*u*2.56E2+(t*t)*1.6E2-(t*t*t)*1.28E2;
  phiu[28] =  s*9.6E1-s*t*2.24E2-s*u*4.48E2+s*(t*t)*1.28E2+(s*s)*t*2.56E2+s*(u*u)*3.84E2+(s*s)*u*5.12E2-
           (s*s)*2.24E2+(s*s*s)*1.28E2+s*t*u*5.12E2;
  phiu[29] =  s*-3.2E1+s*t*3.2E1+s*u*6.4E1-(s*s)*t*1.28E2-(s*s)*u*2.56E2+(s*s)*1.6E2-(s*s*s)*1.28E2;
  phiu[30] =  s*-3.2E1+s*t*3.2E1+s*u*3.2E2-s*(u*u)*3.84E2-(s*s)*u*2.56E2+(s*s)*3.2E1-s*t*u*2.56E2;
  phiu[31] =  s*t*-2.24E2+s*(t*t)*2.56E2+(s*s)*t*2.56E2+s*t*u*2.56E2;
  phiu[32] =  s*t*3.2E1-s*(t*t)*1.28E2;
  phiu[33] =  s*t*3.2E1-(s*s)*t*1.28E2;
  phiu[34] =  s*t*2.56E2-s*(t*t)*2.56E2-(s*s)*t*2.56E2-s*t*u*5.12E2;
}

template<>
void
Lagrange<Simplex,3,4>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  // phiss
  real_t* phiss = phi;
  phiss[0] =  ONETHIRD*1.4E2-t*1.6E2-u*1.6E2-ONETHIRD*s*4.8E2+t*u*2.56E2+ONETHIRD*(s*s)*3.84E2+(t*t)*1.28E2+
           (u*u)*1.28E2+ONETHIRD*s*t*7.68E2+ONETHIRD*s*u*7.68E2;
  phiss[1] =  ONETHIRD*4.4E1-s*9.6E1+ONETHIRD*(s*s)*3.84E2;
  phiss[2] =  0.0;
  phiss[3] =  0.0;
  phiss[4] =  0.0;
  phiss[5] =  0.0;
  phiss[6] =  0.0;
  phiss[7] =  0.0;
  phiss[8] =  u*-3.2E1+(u*u)*1.28E2;
  phiss[9] =  u*-6.4E1+ONETHIRD*s*u*7.68E2;
  phiss[10] =  t*-6.4E1+ONETHIRD*s*t*7.68E2;
  phiss[11] =  t*-3.2E1+(t*t)*1.28E2;
  phiss[12] =  0.0;
  phiss[13] =  0.0;
  phiss[14] =  t*-3.2E1+(t*t)*1.28E2;
  phiss[15] =  t*1.92E2-t*u*2.56E2-(t*t)*2.56E2-ONETHIRD*s*t*7.68E2;
  phiss[16] =  u*1.92E2-t*u*2.56E2-(u*u)*2.56E2-ONETHIRD*s*u*7.68E2;
  phiss[17] =  u*-3.2E1+(u*u)*1.28E2;
  phiss[18] =  0.0;
  phiss[19] =  ONETHIRD*-4.16E2+s*5.76E2+t*3.84E2+u*3.84E2-s*t*7.68E2-s*u*7.68E2-t*u*5.12E2-ONETHIRD*(s*
           s)*1.536E3-(t*t)*2.56E2-(u*u)*2.56E2;
  phiss[20] =  s*-7.68E2-t*2.88E2-u*2.88E2+s*t*7.68E2+s*u*7.68E2+t*u*2.56E2+(s*s)*7.68E2+(t*t)*1.28E2+(u*
           u)*1.28E2+1.52E2;
  phiss[21] =  ONETHIRD*-2.24E2+t*6.4E1+u*6.4E1+ONETHIRD*s*1.344E3-ONETHIRD*(s*s)*1.536E3-ONETHIRD*s*t*
           7.68E2-ONETHIRD*s*u*7.68E2;
  phiss[22] =  t*u*2.56E2;
  phiss[23] =  0.0;
  phiss[24] =  0.0;
  phiss[25] =  t*u*2.56E2;
  phiss[26] =  0.0;
  phiss[27] =  0.0;
  phiss[28] =  u*-4.48E2+s*u*7.68E2+t*u*5.12E2+(u*u)*5.12E2;
  phiss[29] =  u*3.2E2-s*u*7.68E2-t*u*2.56E2-(u*u)*2.56E2;
  phiss[30] =  u*6.4E1-(u*u)*2.56E2;
  phiss[31] =  t*-4.48E2+s*t*7.68E2+t*u*5.12E2+(t*t)*5.12E2;
  phiss[32] =  t*6.4E1-(t*t)*2.56E2;
  phiss[33] =  t*3.2E2-s*t*7.68E2-t*u*2.56E2-(t*t)*2.56E2;
  phiss[34] =  t*u*-5.12E2;

  // phist
  real_t* phist = phi + 35;
  phist[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2-u*1.6E2+s*t*2.56E2+s*u*2.56E2+t*u*2.56E2+ONETHIRD*(s*s)*3.84E2+
           ONETHIRD*(t*t)*3.84E2+(u*u)*1.28E2;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  0.0;
  phist[4] =  0.0;
  phist[5] =  0.0;
  phist[6] =  0.0;
  phist[7] =  0.0;
  phist[8] =  0.0;
  phist[9] =  0.0;
  phist[10] =  ONETHIRD*1.6E1-s*6.4E1+ONETHIRD*(s*s)*3.84E2;
  phist[11] =  s*-3.2E1-t*3.2E1+s*t*2.56E2+4.0;
  phist[12] =  ONETHIRD*1.6E1-t*6.4E1+ONETHIRD*(t*t)*3.84E2;
  phist[13] =  ONETHIRD*-1.6E1+t*6.4E1-ONETHIRD*(t*t)*3.84E2;
  phist[14] =  s*-3.2E1-t*2.88E2-u*3.2E1+s*t*2.56E2+t*u*2.56E2+(t*t)*3.84E2+2.8E1;
  phist[15] =  ONETHIRD*-2.08E2+s*1.92E2+t*3.84E2+u*1.92E2-s*t*5.12E2-s*u*2.56E2-t*u*5.12E2-ONETHIRD*(s*
           s)*3.84E2-(t*t)*3.84E2-(u*u)*1.28E2;
  phist[16] =  u*1.92E2-s*u*2.56E2-t*u*2.56E2-(u*u)*2.56E2;
  phist[17] =  u*-3.2E1+(u*u)*1.28E2;
  phist[18] =  0.0;
  phist[19] =  ONETHIRD*-2.08E2+s*3.84E2+t*1.92E2+u*1.92E2-s*t*5.12E2-s*u*5.12E2-t*u*2.56E2-ONETHIRD*(t*
           t)*3.84E2-(s*s)*3.84E2-(u*u)*1.28E2;
  phist[20] =  s*-2.88E2-t*3.2E1-u*3.2E1+s*t*2.56E2+s*u*2.56E2+(s*s)*3.84E2+2.8E1;
  phist[21] =  ONETHIRD*-1.6E1+s*6.4E1-ONETHIRD*(s*s)*3.84E2;
  phist[22] =  u*-3.2E1+s*u*2.56E2;
  phist[23] =  u*-3.2E1+t*u*2.56E2;
  phist[24] =  u*-3.2E1+(u*u)*1.28E2;
  phist[25] =  u*-2.24E2+s*u*2.56E2+t*u*5.12E2+(u*u)*2.56E2;
  phist[26] =  u*3.2E1-(u*u)*1.28E2;
  phist[27] =  u*3.2E1-t*u*2.56E2;
  phist[28] =  u*-2.24E2+s*u*5.12E2+t*u*2.56E2+(u*u)*2.56E2;
  phist[29] =  u*3.2E1-s*u*2.56E2;
  phist[30] =  u*3.2E1-(u*u)*1.28E2;
  phist[31] =  s*-4.48E2-t*4.48E2-u*2.24E2+s*t*1.024E3+s*u*5.12E2+t*u*5.12E2+(s*s)*3.84E2+(t*t)*3.84E2+
           (u*u)*1.28E2+9.6E1;
  phist[32] =  s*6.4E1+t*3.2E2+u*3.2E1-s*t*5.12E2-t*u*2.56E2-(t*t)*3.84E2-3.2E1;
  phist[33] =  s*3.2E2+t*6.4E1+u*3.2E1-s*t*5.12E2-s*u*2.56E2-(s*s)*3.84E2-3.2E1;
  phist[34] =  u*2.56E2-s*u*5.12E2-t*u*5.12E2-(u*u)*2.56E2;

  // phitt
  real_t* phitt = phi + 70;
  phitt[0] =  ONETHIRD*1.4E2-s*1.6E2-u*1.6E2-ONETHIRD*t*4.8E2+s*u*2.56E2+ONETHIRD*(t*t)*3.84E2+(s*s)*1.28E2+
           (u*u)*1.28E2+ONETHIRD*s*t*7.68E2+ONETHIRD*t*u*7.68E2;
  phitt[1] =  0.0;
  phitt[2] =  ONETHIRD*4.4E1-t*9.6E1+ONETHIRD*(t*t)*3.84E2;
  phitt[3] =  0.0;
  phitt[4] =  u*-6.4E1+ONETHIRD*t*u*7.68E2;
  phitt[5] =  u*-3.2E1+(u*u)*1.28E2;
  phitt[6] =  0.0;
  phitt[7] =  0.0;
  phitt[8] =  0.0;
  phitt[9] =  0.0;
  phitt[10] =  0.0;
  phitt[11] =  s*-3.2E1+(s*s)*1.28E2;
  phitt[12] =  s*-6.4E1+ONETHIRD*s*t*7.68E2;
  phitt[13] =  ONETHIRD*-2.24E2+s*6.4E1+u*6.4E1+ONETHIRD*t*1.344E3-ONETHIRD*(t*t)*1.536E3-ONETHIRD*s*t*
           7.68E2-ONETHIRD*t*u*7.68E2;
  phitt[14] =  s*-2.88E2-t*7.68E2-u*2.88E2+s*t*7.68E2+s*u*2.56E2+t*u*7.68E2+(s*s)*1.28E2+(t*t)*7.68E2+(u*
           u)*1.28E2+1.52E2;
  phitt[15] =  ONETHIRD*-4.16E2+s*3.84E2+t*5.76E2+u*3.84E2-s*t*7.68E2-s*u*5.12E2-t*u*7.68E2-ONETHIRD*(t*
           t)*1.536E3-(s*s)*2.56E2-(u*u)*2.56E2;
  phitt[16] =  u*1.92E2-s*u*2.56E2-(u*u)*2.56E2-ONETHIRD*t*u*7.68E2;
  phitt[17] =  u*-3.2E1+(u*u)*1.28E2;
  phitt[18] =  0.0;
  phitt[19] =  s*1.92E2-s*u*2.56E2-(s*s)*2.56E2-ONETHIRD*s*t*7.68E2;
  phitt[20] =  s*-3.2E1+(s*s)*1.28E2;
  phitt[21] =  0.0;
  phitt[22] =  0.0;
  phitt[23] =  s*u*2.56E2;
  phitt[24] =  0.0;
  phitt[25] =  u*-4.48E2+s*u*5.12E2+t*u*7.68E2+(u*u)*5.12E2;
  phitt[26] =  u*6.4E1-(u*u)*2.56E2;
  phitt[27] =  u*3.2E2-s*u*2.56E2-t*u*7.68E2-(u*u)*2.56E2;
  phitt[28] =  s*u*2.56E2;
  phitt[29] =  0.0;
  phitt[30] =  0.0;
  phitt[31] =  s*-4.48E2+s*t*7.68E2+s*u*5.12E2+(s*s)*5.12E2;
  phitt[32] =  s*3.2E2-s*t*7.68E2-s*u*2.56E2-(s*s)*2.56E2;
  phitt[33] =  s*6.4E1-(s*s)*2.56E2;
  phitt[34] =  s*u*-5.12E2;

  // phisu
  real_t* phisu = phi + 105;
  phisu[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2-u*1.6E2+s*t*2.56E2+s*u*2.56E2+t*u*2.56E2+ONETHIRD*(s*s)*3.84E2+
           ONETHIRD*(u*u)*3.84E2+(t*t)*1.28E2;
  phisu[1] =  0.0;
  phisu[2] =  0.0;
  phisu[3] =  0.0;
  phisu[4] =  0.0;
  phisu[5] =  0.0;
  phisu[6] =  0.0;
  phisu[7] =  ONETHIRD*1.6E1-u*6.4E1+ONETHIRD*(u*u)*3.84E2;
  phisu[8] =  s*-3.2E1-u*3.2E1+s*u*2.56E2+4.0;
  phisu[9] =  ONETHIRD*1.6E1-s*6.4E1+ONETHIRD*(s*s)*3.84E2;
  phisu[10] =  0.0;
  phisu[11] =  0.0;
  phisu[12] =  0.0;
  phisu[13] =  0.0;
  phisu[14] =  t*-3.2E1+(t*t)*1.28E2;
  phisu[15] =  t*1.92E2-s*t*2.56E2-t*u*2.56E2-(t*t)*2.56E2;
  phisu[16] =  ONETHIRD*-2.08E2+s*1.92E2+t*1.92E2+u*3.84E2-s*t*2.56E2-s*u*5.12E2-t*u*5.12E2-ONETHIRD*(s*
           s)*3.84E2-(t*t)*1.28E2-(u*u)*3.84E2;
  phisu[17] =  s*-3.2E1-t*3.2E1-u*2.88E2+s*u*2.56E2+t*u*2.56E2+(u*u)*3.84E2+2.8E1;
  phisu[18] =  ONETHIRD*-1.6E1+u*6.4E1-ONETHIRD*(u*u)*3.84E2;
  phisu[19] =  ONETHIRD*-2.08E2+s*3.84E2+t*1.92E2+u*1.92E2-s*t*5.12E2-s*u*5.12E2-t*u*2.56E2-ONETHIRD*(u*
           u)*3.84E2-(s*s)*3.84E2-(t*t)*1.28E2;
  phisu[20] =  s*-2.88E2-t*3.2E1-u*3.2E1+s*t*2.56E2+s*u*2.56E2+(s*s)*3.84E2+2.8E1;
  phisu[21] =  ONETHIRD*-1.6E1+s*6.4E1-ONETHIRD*(s*s)*3.84E2;
  phisu[22] =  t*-3.2E1+s*t*2.56E2;
  phisu[23] =  t*-3.2E1+(t*t)*1.28E2;
  phisu[24] =  t*-3.2E1+t*u*2.56E2;
  phisu[25] =  t*-2.24E2+s*t*2.56E2+t*u*5.12E2+(t*t)*2.56E2;
  phisu[26] =  t*3.2E1-t*u*2.56E2;
  phisu[27] =  t*3.2E1-(t*t)*1.28E2;
  phisu[28] =  s*-4.48E2-t*2.24E2-u*4.48E2+s*t*5.12E2+s*u*1.024E3+t*u*5.12E2+(s*s)*3.84E2+(t*t)*1.28E2+
           (u*u)*3.84E2+9.6E1;
  phisu[29] =  s*3.2E2+t*3.2E1+u*6.4E1-s*t*2.56E2-s*u*5.12E2-(s*s)*3.84E2-3.2E1;
  phisu[30] =  s*6.4E1+t*3.2E1+u*3.2E2-s*u*5.12E2-t*u*2.56E2-(u*u)*3.84E2-3.2E1;
  phisu[31] =  t*-2.24E2+s*t*5.12E2+t*u*2.56E2+(t*t)*2.56E2;
  phisu[32] =  t*3.2E1-(t*t)*1.28E2;
  phisu[33] =  t*3.2E1-s*t*2.56E2;
  phisu[34] =  t*2.56E2-s*t*5.12E2-t*u*5.12E2-(t*t)*2.56E2;

  // phitu
  real_t* phitu = phi + 140;
  phitu[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2-u*1.6E2+s*t*2.56E2+s*u*2.56E2+t*u*2.56E2+ONETHIRD*(t*t)*3.84E2+
           ONETHIRD*(u*u)*3.84E2+(s*s)*1.28E2;
  phitu[1] =  0.0;
  phitu[2] =  0.0;
  phitu[3] =  0.0;
  phitu[4] =  ONETHIRD*1.6E1-t*6.4E1+ONETHIRD*(t*t)*3.84E2;
  phitu[5] =  t*-3.2E1-u*3.2E1+t*u*2.56E2+4.0;
  phitu[6] =  ONETHIRD*1.6E1-u*6.4E1+ONETHIRD*(u*u)*3.84E2;
  phitu[7] =  0.0;
  phitu[8] =  0.0;
  phitu[9] =  0.0;
  phitu[10] =  0.0;
  phitu[11] =  0.0;
  phitu[12] =  0.0;
  phitu[13] =  ONETHIRD*-1.6E1+t*6.4E1-ONETHIRD*(t*t)*3.84E2;
  phitu[14] =  s*-3.2E1-t*2.88E2-u*3.2E1+s*t*2.56E2+t*u*2.56E2+(t*t)*3.84E2+2.8E1;
  phitu[15] =  ONETHIRD*-2.08E2+s*1.92E2+t*3.84E2+u*1.92E2-s*t*5.12E2-s*u*2.56E2-t*u*5.12E2-ONETHIRD*(u*
           u)*3.84E2-(s*s)*1.28E2-(t*t)*3.84E2;
  phitu[16] =  ONETHIRD*-2.08E2+s*1.92E2+t*1.92E2+u*3.84E2-s*t*2.56E2-s*u*5.12E2-t*u*5.12E2-ONETHIRD*(t*
           t)*3.84E2-(s*s)*1.28E2-(u*u)*3.84E2;
  phitu[17] =  s*-3.2E1-t*3.2E1-u*2.88E2+s*u*2.56E2+t*u*2.56E2+(u*u)*3.84E2+2.8E1;
  phitu[18] =  ONETHIRD*-1.6E1+u*6.4E1-ONETHIRD*(u*u)*3.84E2;
  phitu[19] =  s*1.92E2-s*t*2.56E2-s*u*2.56E2-(s*s)*2.56E2;
  phitu[20] =  s*-3.2E1+(s*s)*1.28E2;
  phitu[21] =  0.0;
  phitu[22] =  s*-3.2E1+(s*s)*1.28E2;
  phitu[23] =  s*-3.2E1+s*t*2.56E2;
  phitu[24] =  s*-3.2E1+s*u*2.56E2;
  phitu[25] =  s*-2.24E2-t*4.48E2-u*4.48E2+s*t*5.12E2+s*u*5.12E2+t*u*1.024E3+(s*s)*1.28E2+(t*t)*3.84E2+
           (u*u)*3.84E2+9.6E1;
  phitu[26] =  s*3.2E1+t*6.4E1+u*3.2E2-s*u*2.56E2-t*u*5.12E2-(u*u)*3.84E2-3.2E1;
  phitu[27] =  s*3.2E1+t*3.2E2+u*6.4E1-s*t*2.56E2-t*u*5.12E2-(t*t)*3.84E2-3.2E1;
  phitu[28] =  s*-2.24E2+s*t*2.56E2+s*u*5.12E2+(s*s)*2.56E2;
  phitu[29] =  s*3.2E1-(s*s)*1.28E2;
  phitu[30] =  s*3.2E1-s*u*2.56E2;
  phitu[31] =  s*-2.24E2+s*t*5.12E2+s*u*2.56E2+(s*s)*2.56E2;
  phitu[32] =  s*3.2E1-s*t*2.56E2;
  phitu[33] =  s*3.2E1-(s*s)*1.28E2;
  phitu[34] =  s*2.56E2-s*t*5.12E2-s*u*5.12E2-(s*s)*2.56E2;

  // phiuu
  real_t* phiuu = phi + 175;
  phiuu[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2-ONETHIRD*u*4.8E2+s*t*2.56E2+ONETHIRD*(u*u)*3.84E2+(s*s)*1.28E2+
           (t*t)*1.28E2+ONETHIRD*s*u*7.68E2+ONETHIRD*t*u*7.68E2;
  phiuu[1] =  0.0;
  phiuu[2] =  0.0;
  phiuu[3] =  ONETHIRD*4.4E1-u*9.6E1+ONETHIRD*(u*u)*3.84E2;
  phiuu[4] =  0.0;
  phiuu[5] =  t*-3.2E1+(t*t)*1.28E2;
  phiuu[6] =  t*-6.4E1+ONETHIRD*t*u*7.68E2;
  phiuu[7] =  s*-6.4E1+ONETHIRD*s*u*7.68E2;
  phiuu[8] =  s*-3.2E1+(s*s)*1.28E2;
  phiuu[9] =  0.0;
  phiuu[10] =  0.0;
  phiuu[11] =  0.0;
  phiuu[12] =  0.0;
  phiuu[13] =  0.0;
  phiuu[14] =  t*-3.2E1+(t*t)*1.28E2;
  phiuu[15] =  t*1.92E2-s*t*2.56E2-(t*t)*2.56E2-ONETHIRD*t*u*7.68E2;
  phiuu[16] =  ONETHIRD*-4.16E2+s*3.84E2+t*3.84E2+u*5.76E2-s*t*5.12E2-s*u*7.68E2-t*u*7.68E2-ONETHIRD*(u*
           u)*1.536E3-(s*s)*2.56E2-(t*t)*2.56E2;
  phiuu[17] =  s*-2.88E2-t*2.88E2-u*7.68E2+s*t*2.56E2+s*u*7.68E2+t*u*7.68E2+(s*s)*1.28E2+(t*t)*1.28E2+(u*
           u)*7.68E2+1.52E2;
  phiuu[18] =  ONETHIRD*-2.24E2+s*6.4E1+t*6.4E1+ONETHIRD*u*1.344E3-ONETHIRD*(u*u)*1.536E3-ONETHIRD*s*u*
           7.68E2-ONETHIRD*t*u*7.68E2;
  phiuu[19] =  s*1.92E2-s*t*2.56E2-(s*s)*2.56E2-ONETHIRD*s*u*7.68E2;
  phiuu[20] =  s*-3.2E1+(s*s)*1.28E2;
  phiuu[21] =  0.0;
  phiuu[22] =  0.0;
  phiuu[23] =  0.0;
  phiuu[24] =  s*t*2.56E2;
  phiuu[25] =  t*-4.48E2+s*t*5.12E2+t*u*7.68E2+(t*t)*5.12E2;
  phiuu[26] =  t*3.2E2-s*t*2.56E2-t*u*7.68E2-(t*t)*2.56E2;
  phiuu[27] =  t*6.4E1-(t*t)*2.56E2;
  phiuu[28] =  s*-4.48E2+s*t*5.12E2+s*u*7.68E2+(s*s)*5.12E2;
  phiuu[29] =  s*6.4E1-(s*s)*2.56E2;
  phiuu[30] =  s*3.2E2-s*t*2.56E2-s*u*7.68E2-(s*s)*2.56E2;
  phiuu[31] =  s*t*2.56E2;
  phiuu[32] =  0.0;
  phiuu[33] =  0.0;
  phiuu[34] =  s*t*-5.12E2;
}

/*
 * pentatopes, p = 5
*/
template<>
void
Lagrange<Simplex,3,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,3,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,3,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
