#include "common/error.h"
#include "element/basis.h"

namespace avro
{

/*
 * tetrahedron, p = 1
*/

template<>
void
Bernstein<Simplex,3,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t b1 = 1 - s - t - u;
  real_t b2 = s;
  real_t b3 = t;
  real_t b4 = u;

  phi[0] = b1;
  phi[1] = b2;
  phi[2] = b3;
  phi[3] = b4;
}

template<>
void
Bernstein<Simplex,3,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = -1;
  phis[1] = 1;
  phis[2] = 0;
  phis[3] = 0;

  real_t* phit = phi + 4;
  phit[0] = -1;
  phit[1] = 0;
  phit[2] = 1;
  phit[3] = 0;

  real_t* phiu = phi + 8;
  phiu[0] = -1;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = 1;
}

template<>
void
Bernstein<Simplex,3,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 24; i++)
    phi[i] = 0;
}

/*
 * tetrahedron, p = 2
*/
template<>
void
Bernstein<Simplex,3,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t b1 = 1 - s - t - u;
  real_t b2 = s;
  real_t b3 = t;
  real_t b4 = u;

  phi[0] = b1*b1;
  phi[1] = b2*b2;
  phi[2] = b3*b3;
  phi[3] = b4*b4;
  phi[4] = 2*b3*b4;
  phi[5] = 2*b2*b4;
  phi[6] = 2*b2*b3;
  phi[7] = 2*b1*b3;
  phi[8] = 2*b1*b4;
  phi[9] = 2*b1*b2;
}

template<>
void
Bernstein<Simplex,3,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t* phis = phi;
  phis[0] = -2*(1 - s - t - u);
  phis[1] = 2*s;
  phis[2] = 0;
  phis[3] = 0;
  phis[4] = 0;
  phis[5] = 2*u;
  phis[6] = 2*t;
  phis[7] = -2*t;
  phis[8] = -2*u;
  phis[9] = -2* s +2 *(1 - s - t - u);

  real_t* phit = phi + 10;
  phit[0] = -2*(1-s-t-u);
  phit[1] = 0;
  phit[2] = 2*t;
  phit[3] = 0;
  phit[4] = 2*u;
  phit[5] = 0;
  phit[6] = 2*s;
  phit[7] = -2*t+2*(1-s-t-u);
  phit[8] = -2*u;
  phit[9] = -2*s;

  real_t* phiu = phi + 20;
  phiu[0] = -2*(1-s-t-u);
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = 2*u;
  phiu[4] = 2*t;
  phiu[5] = 2*s;
  phiu[6] = 0;
  phiu[7] = -2*t;
  phiu[8] = 2*(1-s-t-u)-2*u;
  phiu[9] = -2*s;

}

template<>
void
Bernstein<Simplex,3,2>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * tetrahedron, p = 3
*/
template<>
void
Bernstein<Simplex,3,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t b1 = 1 - s - t - u;
  real_t b2 = s;
  real_t b3 = t;
  real_t b4 = u;

  phi[0]  =  b1*b1*b1;
  phi[1]  =  b2*b2*b2;
  phi[2]  =  b3*b3*b3;
  phi[3]  =  b4*b4*b4;
  phi[4]  = 3*b3*b3*b4;
  phi[5]  = 3*b3*b4*b4;
  phi[6]  = 3*b2*b4*b4;
  phi[7]  = 3*b2*b2*b4;
  phi[8]  = 3*b2*b2*b3;
  phi[9]  = 3*b2*b3*b3;
  phi[10] = 3*b1*b3*b3;
  phi[11] = 3*b1*b1*b3;
  phi[12] = 3*b1*b1*b4;
  phi[13] = 3*b1*b4*b4;
  phi[14] = 3*b1*b1*b2;
  phi[15] = 3*b1*b2*b2;
  phi[16] = 6*b2*b3*b4;
  phi[17] = 6*b1*b3*b4;
  phi[18] = 6*b1*b2*b4;
  phi[19] = 6*b1*b2*b3;
}

template<>
void
Bernstein<Simplex,3,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t* phis = phi;
  phis[0]  = -3*(1-s-t-u)*(1-s-t-u);
  phis[1]  = 3*s*s;
  phis[2]  = 0;
  phis[3]  = 0;
  phis[4]  = 0;
  phis[5]  = 0;
  phis[6]  = 3*u*u;
  phis[7]  = 6*s*u;
  phis[8]  = 6*s*t;
  phis[9]  = 3*t*t;
  phis[10] = -3*t*t;
  phis[11] = -6*t*(1-s-t-u);
  phis[12] = -6*(1-s-t-u)*u;
  phis[13] = -3*u*u;
  phis[14] = -6*s*(1-s-t-u)+3*(1-s-t-u)*(1-s-t-u);
  phis[15] = -3*s*s+6*s*(1-s-t-u);
  phis[16] = 6*t*u;
  phis[17] = -6*t*u;
  phis[18] = -6*s*u+6*(1-s-t-u)*u;
  phis[19] = -6*s*t+6*t*(1-s-t-u);

  real_t* phit = phi + 20;
  phit[0]  = -3*(1-s-t-u)*(1-s-t-u);
  phit[1]  = 0;
  phit[2]  = 3*t*t;
  phit[3]  = 0;
  phit[4]  = 6*t*u;
  phit[5]  = 3*u*u;
  phit[6]  = 0;
  phit[7]  = 0;
  phit[8]  = 3*s*s;
  phit[9]  = 6*s*t;
  phit[10] = -3*t*t+6*t*(1-s-t-u);
  phit[11] = -6*t*(1-s-t-u)+3*(1-s-t-u)*(1-s-t-u);
  phit[12] = -6*(1-s-t-u)*u;
  phit[13] = -3*u*u;
  phit[14] = -6*s*(1-s-t-u);
  phit[15] = -3*s*s;
  phit[16] = 6*s*u;
  phit[17] = -6*t*u+6*(1-s-t-u)*u;
  phit[18] = -6*s*u;
  phit[19] = -6*s*t+6*s*(1-s-t-u);

  real_t* phiu = phi + 40;
  phiu[0] = -3*(1-s-t-u)*(1-s-t-u);
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = 3*u*u;
  phiu[4] = 3*t*t;
  phiu[5] = 6*t*u;
  phiu[6] = 6*s*u;
  phiu[6] = 3*s*s;
  phiu[8] = 0;
  phiu[9] = 0;
  phiu[10] = -3*t*t;
  phiu[11] = -6*t*(1-s-t-u);
  phiu[12] = 3*(1-s-t-u)*(1-s-t-u)-6*(1-s-t-u)*u;
  phiu[13] = 6*(1-s-t-u)*u-3*u*u;
  phiu[14] = -6*s*(1-s-t-u);
  phiu[15] = -3*s*s;
  phiu[16] = 6*s*t;
  phiu[17] = 6*t*(1-s-t-u)-6*t*u;
  phiu[18] = 6*s*(1-s-t-u)-6*s*u;
  phiu[19] = -6*s*t;
}

template<>
void
Bernstein<Simplex,3,3>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * tetrahedron, p = 4
*/
template<>
void
Bernstein<Simplex,3,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t b1 = 1 - s - t - u;
  real_t b2 = s;
  real_t b3 = t;
  real_t b4 = u;

  phi[0] = b1*b1*b1*b1;
  phi[1] = b2*b2*b2*b2;
  phi[2] = b3*b3*b3*b3;
  phi[3] = b4*b4*b4*b4;
  phi[4] = 4*b1*b1*b1*b2;
  phi[5] = 6*b1*b1*b2*b2;
  phi[6] = 4*b1*b2*b2*b2;
  phi[7] = 4*b2*b2*b2*b3;
  phi[8] = 6*b2*b2*b3*b3;
  phi[9] = 4*b2*b3*b3*b3;
  phi[10] = 4*b3*b3*b3*b1;
  phi[11] = 6*b3*b3*b1*b1;
  phi[12] = 4*b3*b1*b1*b1;
  phi[13] = 4*b3*b3*b3*b4;
  phi[14] = 6*b3*b3*b4*b4;
  phi[15] = 4*b3*b4*b4*b4;
  phi[16] = 4*b1*b1*b1*b4;
  phi[17] = 6*b1*b1*b4*b4;
  phi[18] = 4*b1*b4*b4*b4;
  phi[19] = 4*b2*b2*b2*b4;
  phi[20] = 6*b2*b2*b4*b4;
  phi[21] = 4*b2*b4*b4*b4;
  phi[22] = 12*b1*b1*b2*b3;
  phi[23] = 12*b1*b2*b2*b3;
  phi[24] = 12*b1*b2*b3*b3;
  phi[25] = 12*b1*b1*b3*b4;
  phi[26] = 12*b1*b4*b4*b3;
  phi[27] = 12*b1*b4*b3*b3;
  phi[28] = 12*b1*b1*b2*b4;
  phi[29] = 12*b1*b2*b2*b4;
  phi[30] = 12*b1*b2*b4*b4;
  phi[31] = 12*b2*b2*b3*b4;
  phi[32] = 12*b2*b3*b3*b4;
  phi[33] = 12*b2*b3*b4*b4;
  phi[34] = 24*b1*b2*b3*b4;
}

template<>
void
Bernstein<Simplex,3,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  // phis
  real_t* phis = phi;
  phis[0] = -4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u);
  phis[1] = 4 * s * s * s;
  phis[2] = 0;
  phis[3] = 0;
  phis[4] = -12 * s * (1 - s - t - u) * (1 - s - t - u) + 4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u);
  phis[5] = -12 * s * s * (1 - s - t - u) + 12 * s * (1 - s - t - u) * (1 - s - t - u);
  phis[6] = -4 * s * s * s + 12 * s * s * (1 - s - t - u);
  phis[7] = 12 * s * s * t;
  phis[8] = 12 * s * t * t;
  phis[9] = 4 * t * t * t;
  phis[10] = -4 * t * t * t;
  phis[11] = -12 * t * t * (1 - s - t - u);
  phis[12] = -12 * t * (1 - s - t - u) * (1 - s - t - u);
  phis[13] = 0;
  phis[14] = 0;
  phis[15] = 0;
  phis[16] = -12 * (1 - s - t - u) * (1 - s - t - u) * u;
  phis[17] = -12 * (1 - s - t - u) * u * u;
  phis[18] = -4 * u * u * u;
  phis[19] = 12 * s * s * u;
  phis[20] = 12 * s * u * u;
  phis[21] = 4 * u * u * u;
  phis[22] = -24 * s * t * (1 - s - t - u) + 12 * t * (1 - s - t - u) * (1 - s - t - u);
  phis[23] = -12 * s * s * t + 24 * s * t * (1 - s - t - u);
  phis[24] = -12 * s * t * t + 12 * t * t * (1 - s - t - u);
  phis[25] = -24 * t * (1 - s - t - u) * u;
  phis[26] = -12 * t * u * u;
  phis[27] = -12 * t * t * u;
  phis[28] = -24 * s * (1 - s - t - u) * u + 12 * (1 - s - t - u) * (1 - s - t - u) * u;
  phis[29] = -12 * s * s * u + 24 * s * (1 - s - t - u) * u;
  phis[30] = -12 * s * u * u + 12 * (1 - s - t - u) * u * u;
  phis[31] = 24 * s * t * u;
  phis[32] = 12 * t * t * u;
  phis[33] = 12 * t * u * u;
  phis[34] = -24 * s * t * u + 24 * t * (1 - s - t - u) * u;

  // phit
  real_t *phit = phi + 35;
  phit[0] = -4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u);
  phit[1] = 0;
  phit[2] = 4 * t * t * t;
  phit[3] = 0;
  phit[4] = -12 * s * (1 - s - t - u) * (1 - s - t - u);
  phit[5] = -12 * s * s * (1 - s - t - u);
  phit[6] = -4 * s * s * s;
  phit[7] = 4 * s * s * s;
  phit[8] = 12 * s * s * t;
  phit[9] = 12 * s * t * t;
  phit[10] = -4 * t * t * t + 12 * t * t * (1 - s - t - u);
  phit[11] = -12 * t * t * (1 - s - t - u) + 12 * t * (1 - s - t - u) * (1 - s - t - u);
  phit[12] = -12 * t * (1 - s - t - u) * (1 - s - t - u) + 4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u);
  phit[13] = 12 * t * t * u;
  phit[14] = 12 * t * u * u;
  phit[15] = 4 * u * u * u;
  phit[16] = -12 * (1 - s - t - u) * (1 - s - t - u) * u;
  phit[17] = -12 * (1 - s - t - u) * u * u;
  phit[18] = -4 * u * u * u;
  phit[19] = 0;
  phit[20] = 0;
  phit[21] = 0;
  phit[22] = -24 * s * t * (1 - s - t - u) + 12 * s * (1 - s - t - u) * (1 - s - t - u);
  phit[23] = -12 * s * s * t + 12 * s * s * (1 - s - t - u);
  phit[24] = -12 * s * t * t + 24 * s * t * (1 - s - t - u);
  phit[25] = -24 * t * (1 - s - t - u) * u + 12 * (1 - s - t - u) * (1 - s - t - u) * u;
  phit[26] = -12 * t * u * u + 12 * (1 - s - t - u) * u * u;
  phit[27] = -12 * t * t * u + 24 * t * (1 - s - t - u) * u;
  phit[28] = -24 * s * (1 - s - t - u) * u;
  phit[29] = -12 * s * s * u;
  phit[30] = -12 * s * u * u;
  phit[31] = 12 * s * s * u;
  phit[32] = 24 * s * t * u;
  phit[33] = 12 * s * u * u;
  phit[34] = -24 * s * t * u + 24 * s * (1 - s - t - u) * u;

  // phiu
  real_t* phiu = phi + 70;
  phiu[0] = -4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u);
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = 4 * u * u * u;
  phiu[4] = -12 * s * (1 - s - t - u) * (1 - s - t - u);
  phiu[5] = -12 * s * s * (1 - s - t - u);
  phiu[6] = -4 * s * s * s;
  phiu[7] = 0;
  phiu[8] = 0;
  phiu[9] = 0;
  phiu[10] = -4 * t * t * t;
  phiu[11] = -12 * t * t * (1 - s - t - u);
  phiu[12] = -12 * t * (1 - s - t - u) * (1 - s - t - u);
  phiu[13] = 4 * t * t * t;
  phiu[14] = 12 * t * t * u;
  phiu[15] = 12 * t * u * u;
  phiu[16] = 4 * (1 - s - t - u) * (1 - s - t - u) * (1 - s - t - u) - 12 * (1 - s - t - u) * (1 - s - t - u) * u;
  phiu[17] = 12 * (1 - s - t - u) * (1 - s - t - u) * u - 12 * (1 - s - t - u) * u * u;
  phiu[18] = 12 * (1 - s - t - u) * u * u - 4 * u * u * u;
  phiu[19] = 4 * s * s * s;
  phiu[20] = 12 * s * s * u;
  phiu[21] = 12 * s * u * u;
  phiu[22] = -24 * s * t * (1 - s - t - u);
  phiu[23] = -12 * s * s * t;
  phiu[24] = -12 * s * t * t;
  phiu[25] = 12 * t * (1 - s - t - u) * (1 - s - t - u) - 24 * t * (1 - s - t - u) * u;
  phiu[26] = 24 * t * (1 - s - t - u) * u - 12 * t * u * u;
  phiu[27] = 12 * t * t * (1 - s - t - u) - 12 * t * t * u;
  phiu[28] = 12 * s * (1 - s - t - u) * (1 - s - t - u) - 24 * s * (1 - s - t - u) * u;
  phiu[29] = 12 * s * s * (1 - s - t - u) - 12 * s * s * u;
  phiu[30] = 24 * s * (1 - s - t - u) * u - 12 * s * u * u;
  phiu[31] = 12 * s * s * t;
  phiu[32] = 12 * s * t * t;
  phiu[33] = 24 * s * t * u;
  phiu[34] = 24 * s * t * (1 - s - t - u) - 24 * s * t * u;
}

template<>
void
Bernstein<Simplex,3,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * tetrahedron, p = 5
*/
template<>
void
Bernstein<Simplex,3,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,3,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,3,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
