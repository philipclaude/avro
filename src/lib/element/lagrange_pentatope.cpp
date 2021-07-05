#include "common/error.h"
#include "element/basis.h"

namespace avro
{

/*
 * pentatopes, p = 1
*/
template<>
void
Lagrange<Simplex,4,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  phi[0]= 1 - s - t - u - v; // 1 at node 0 (s= t= u= v= 0)
  phi[1]=     s            ; // 1 at node 1 (s= 1, t= u= v= 0)
  phi[2]=         t        ; // 1 at node 2 (t= 1, s= u= v= 0)
  phi[3]=             u    ; // 1 at node 3 (u= 1, s= t= v= 0)
  phi[4]=                 v; // 1 at node 4 (v= 1, s= t= u= 0)
}

template<>
void
Lagrange<Simplex,4,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0]= -1;
  phis[1]= 1;
  phis[2]= 0;
  phis[3]= 0;
  phis[4]= 0;

  real_t* phit = phi + 5;
  phit[0]= -1;
  phit[1]= 0;
  phit[2]= 1;
  phit[3]= 0;
  phit[4]= 0;

  real_t* phiu = phi + 10;
  phiu[0]= -1;
  phiu[1]= 0;
  phiu[2]= 0;
  phiu[3]= 1;
  phiu[4]= 0;

  real_t *phiv = phi + 15;
  phiv[0]= -1;
  phiv[1]= 0;
  phiv[2]= 0;
  phiv[3]= 0;
  phiv[4]= 1;
}

template<>
void
Lagrange<Simplex,4,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 50; i++)
    phi[i] = 0;
}

/*
 * pentatopes, p = 2
*/
template<>
void
Lagrange<Simplex,4,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  phi[0]= s*-3.0-t*3.0-u*3.0-v*3.0+s*t*4.0+s*u*4.0+s*v*4.0+t*u*4.0+t*v*4.0+u*v*4.0+(s*s)*2.0+(t*t)*2.0+
      (u*u)*2.0+(v*v)*2.0+1.0;
  phi[1]= s*(s*2.0-1.0);
  phi[2]= t*(t*2.0-1.0);
  phi[3]= u*(u*2.0-1.0);
  phi[4]= v*(v*2.0-1.0);
  phi[5]= s*(s+t+u+v-1.0)*-4.0;
  phi[6]= t*(s+t+u+v-1.0)*-4.0;
  phi[7]= u*(s+t+u+v-1.0)*-4.0;
  phi[8]= v*(s+t+u+v-1.0)*-4.0;
  phi[9]= s*t*4.0;
  phi[10]= s*u*4.0;
  phi[11]= s*v*4.0;
  phi[12]= t*u*4.0;
  phi[13]= t*v*4.0;
  phi[14]= u*v*4.0;
}

template<>
void
Lagrange<Simplex,4,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  real_t* phis = phi;
  phis[0]= s*4.0+t*4.0+u*4.0+v*4.0-3.0;
  phis[1]= s*4.0-1.0;
  phis[2]= 0.0;
  phis[3]= 0.0;
  phis[4]= 0.0;
  phis[5]= s*-8.0-t*4.0-u*4.0-v*4.0+4.0;
  phis[6]= t*-4.0;
  phis[7]= u*-4.0;
  phis[8]= v*-4.0;
  phis[9]= t*4.0;
  phis[10]= u*4.0;
  phis[11]= v*4.0;
  phis[12]= 0.0;
  phis[13]= 0.0;
  phis[14]= 0.0;

  real_t* phit = phi + 15;
  phit[0]= s*4.0+t*4.0+u*4.0+v*4.0-3.0;
  phit[1]= 0.0;
  phit[2]= t*4.0-1.0;
  phit[3]= 0.0;
  phit[4]= 0.0;
  phit[5]= s*-4.0;
  phit[6]= s*-4.0-t*8.0-u*4.0-v*4.0+4.0;
  phit[7]= u*-4.0;
  phit[8]= v*-4.0;
  phit[9]= s*4.0;
  phit[10]= 0.0;
  phit[11]= 0.0;
  phit[12]= u*4.0;
  phit[13]= v*4.0;
  phit[14]= 0.0;

  real_t* phiu = phi + 30;
  phiu[0]= s*4.0+t*4.0+u*4.0+v*4.0-3.0;
  phiu[1]= 0.0;
  phiu[2]= 0.0;
  phiu[3]= u*4.0-1.0;
  phiu[4]= 0.0;
  phiu[5]= s*-4.0;
  phiu[6]= t*-4.0;
  phiu[7]= s*-4.0-t*4.0-u*8.0-v*4.0+4.0;
  phiu[8]= v*-4.0;
  phiu[9]= 0.0;
  phiu[10]= s*4.0;
  phiu[11]= 0.0;
  phiu[12]= t*4.0;
  phiu[13]= 0.0;
  phiu[14]= v*4.0;


  real_t *phiv = phi + 45;
  phiv[0]= s*4.0+t*4.0+u*4.0+v*4.0-3.0;
  phiv[1]= 0.0;
  phiv[2]= 0.0;
  phiv[3]= 0.0;
  phiv[4]= v*4.0-1.0;
  phiv[5]= s*-4.0;
  phiv[6]= t*-4.0;
  phiv[7]= u*-4.0;
  phiv[8]= s*-4.0-t*4.0-u*4.0-v*8.0+4.0;
  phiv[9]= 0.0;
  phiv[10]= 0.0;
  phiv[11]= s*4.0;
  phiv[12]= 0.0;
  phiv[13]= t*4.0;
  phiv[14]= u*4.0;
}

template<>
void
Lagrange<Simplex,4,2>::hess( const real_t* x , real_t* phi ) {

  // phiss
  real_t *phiss = phi;
  phiss[0]= 4.0;
  phiss[1]= 4.0;
  phiss[2]= 0.0;
  phiss[3]= 0.0;
  phiss[4]= 0.0;
  phiss[5]= -8.0;
  phiss[6]= 0.0;
  phiss[7]= 0.0;
  phiss[8]= 0.0;
  phiss[9]= 0.0;
  phiss[10]= 0.0;
  phiss[11]= 0.0;
  phiss[12]= 0.0;
  phiss[13]= 0.0;
  phiss[14]= 0.0;

  // phist
  real_t *phist = phi + 15;
  phist[0]= 4.0;
  phist[1]= 0.0;
  phist[2]= 0.0;
  phist[3]= 0.0;
  phist[4]= 0.0;
  phist[5]= -4.0;
  phist[6]= -4.0;
  phist[7]= 0.0;
  phist[8]= 0.0;
  phist[9]= 4.0;
  phist[10]= 0.0;
  phist[11]= 0.0;
  phist[12]= 0.0;
  phist[13]= 0.0;
  phist[14]= 0.0;

  // phitt
  real_t *phitt = phi + 30;
  phitt[0]= 4.0;
  phitt[1]= 0.0;
  phitt[2]= 4.0;
  phitt[3]= 0.0;
  phitt[4]= 0.0;
  phitt[5]= 0.0;
  phitt[6]= -8.0;
  phitt[7]= 0.0;
  phitt[8]= 0.0;
  phitt[9]= 0.0;
  phitt[10]= 0.0;
  phitt[11]= 0.0;
  phitt[12]= 0.0;
  phitt[13]= 0.0;
  phitt[14]= 0.0;

  // phisu
  real_t* phisu = phi + 45;
  phisu[0]= 4.0;
  phisu[1]= 0.0;
  phisu[2]= 0.0;
  phisu[3]= 0.0;
  phisu[4]= 0.0;
  phisu[5]= -4.0;
  phisu[6]= 0.0;
  phisu[7]= -4.0;
  phisu[8]= 0.0;
  phisu[9]= 0.0;
  phisu[10]= 4.0;
  phisu[11]= 0.0;
  phisu[12]= 0.0;
  phisu[13]= 0.0;
  phisu[14]= 0.0;

  // phitu
  real_t *phitu = phi + 60;
  phitu[0]= 4.0;
  phitu[1]= 0.0;
  phitu[2]= 0.0;
  phitu[3]= 0.0;
  phitu[4]= 0.0;
  phitu[5]= 0.0;
  phitu[6]= -4.0;
  phitu[7]= -4.0;
  phitu[8]= 0.0;
  phitu[9]= 0.0;
  phitu[10]= 0.0;
  phitu[11]= 0.0;
  phitu[12]= 4.0;
  phitu[13]= 0.0;
  phitu[14]= 0.0;

  // phiuu
  real_t *phiuu = phi + 75;
  phiuu[0]= 4.0;
  phiuu[1]= 0.0;
  phiuu[2]= 0.0;
  phiuu[3]= 4.0;
  phiuu[4]= 0.0;
  phiuu[5]= 0.0;
  phiuu[6]= 0.0;
  phiuu[7]= -8.0;
  phiuu[8]= 0.0;
  phiuu[9]= 0.0;
  phiuu[10]= 0.0;
  phiuu[11]= 0.0;
  phiuu[12]= 0.0;
  phiuu[13]= 0.0;
  phiuu[14]= 0.0;

  // phisv
  real_t *phisv = phi + 90;
  phisv[0]= 4.0;
  phisv[1]= 0.0;
  phisv[2]= 0.0;
  phisv[3]= 0.0;
  phisv[4]= 0.0;
  phisv[5]= -4.0;
  phisv[6]= 0.0;
  phisv[7]= 0.0;
  phisv[8]= -4.0;
  phisv[9]= 0.0;
  phisv[10]= 0.0;
  phisv[11]= 4.0;
  phisv[12]= 0.0;
  phisv[13]= 0.0;
  phisv[14]= 0.0;

  // phitv
  real_t *phitv = phi + 105;
  phitv[0]= 4.0;
  phitv[1]= 0.0;
  phitv[2]= 0.0;
  phitv[3]= 0.0;
  phitv[4]= 0.0;
  phitv[5]= 0.0;
  phitv[6]= -4.0;
  phitv[7]= 0.0;
  phitv[8]= -4.0;
  phitv[9]= 0.0;
  phitv[10]= 0.0;
  phitv[11]= 0.0;
  phitv[12]= 0.0;
  phitv[13]= 4.0;
  phitv[14]= 0.0;

  // phiuv
  real_t *phiuv = phi + 120;
  phiuv[0]= 4.0;
  phiuv[1]= 0.0;
  phiuv[2]= 0.0;
  phiuv[3]= 0.0;
  phiuv[4]= 0.0;
  phiuv[5]= 0.0;
  phiuv[6]= 0.0;
  phiuv[7]= -4.0;
  phiuv[8]= -4.0;
  phiuv[9]= 0.0;
  phiuv[10]= 0.0;
  phiuv[11]= 0.0;
  phiuv[12]= 0.0;
  phiuv[13]= 0.0;
  phiuv[14]= 4.0;

  // phivv
  real_t *phivv = phi + 135;
  phivv[0]= 4.0;
  phivv[1]= 0.0;
  phivv[2]= 0.0;
  phivv[3]= 0.0;
  phivv[4]= 4.0;
  phivv[5]= 0.0;
  phivv[6]= 0.0;
  phivv[7]= 0.0;
  phivv[8]= -8.0;
  phivv[9]= 0.0;
  phivv[10]= 0.0;
  phivv[11]= 0.0;
  phivv[12]= 0.0;
  phivv[13]= 0.0;
  phivv[14]= 0.0;
}

/*
 * pentatopes, p = 3
*/
template<>
void
Lagrange<Simplex,4,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  phi[0]= s*(-1.1E1/2.0)-t*(1.1E1/2.0)-u*(1.1E1/2.0)-v*(1.1E1/2.0)+s*t*1.8E1+s*u*1.8E1+s*v*1.8E1+t*u*
             1.8E1+t*v*1.8E1+u*v*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)-s*(u*u)*(2.7E1/2.0)-(s*s)*u*(2.7E1/
             2.0)-s*(v*v)*(2.7E1/2.0)-t*(u*u)*(2.7E1/2.0)-(s*s)*v*(2.7E1/2.0)-(t*t)*u*(2.7E1/2.0)-t*(v*v)*(2.7E1/2.0)-
             (t*t)*v*(2.7E1/2.0)-u*(v*v)*(2.7E1/2.0)-(u*u)*v*(2.7E1/2.0)+(s*s)*9.0-(s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*
             t*t)*(9.0/2.0)+(u*u)*9.0-(u*u*u)*(9.0/2.0)+(v*v)*9.0-(v*v*v)*(9.0/2.0)-s*t*u*2.7E1-s*t*v*2.7E1-s*u*v*
             2.7E1-t*u*v*2.7E1+1.0;
  phi[1]= (s*(s*-9.0+(s*s)*9.0+2.0))/2.0;
  phi[2]= (t*(t*-9.0+(t*t)*9.0+2.0))/2.0;
  phi[3]= (u*(u*-9.0+(u*u)*9.0+2.0))/2.0;
  phi[4]= (v*(v*-9.0+(v*v)*9.0+2.0))/2.0;
  phi[5]= s*(s*-5.0-t*5.0-u*5.0-v*5.0+s*t*6.0+s*u*6.0+s*v*6.0+t*u*6.0+t*v*6.0+u*v*6.0+(s*s)*3.0+(t*t)*
             3.0+(u*u)*3.0+(v*v)*3.0+2.0)*(9.0/2.0);
  phi[6]= s*(s*3.0-1.0)*(s+t+u+v-1.0)*(-9.0/2.0);
  phi[7]= t*(s*-5.0-t*5.0-u*5.0-v*5.0+s*t*6.0+s*u*6.0+s*v*6.0+t*u*6.0+t*v*6.0+u*v*6.0+(s*s)*3.0+(t*t)*
             3.0+(u*u)*3.0+(v*v)*3.0+2.0)*(9.0/2.0);
  phi[8]= t*(t*3.0-1.0)*(s+t+u+v-1.0)*(-9.0/2.0);
  phi[9]= u*(s*-5.0-t*5.0-u*5.0-v*5.0+s*t*6.0+s*u*6.0+s*v*6.0+t*u*6.0+t*v*6.0+u*v*6.0+(s*s)*3.0+(t*t)*
             3.0+(u*u)*3.0+(v*v)*3.0+2.0)*(9.0/2.0);
  phi[10]= u*(u*3.0-1.0)*(s+t+u+v-1.0)*(-9.0/2.0);
  phi[11]= v*(s*-5.0-t*5.0-u*5.0-v*5.0+s*t*6.0+s*u*6.0+s*v*6.0+t*u*6.0+t*v*6.0+u*v*6.0+(s*s)*3.0+(t*t)*
             3.0+(u*u)*3.0+(v*v)*3.0+2.0)*(9.0/2.0);
  phi[12]= v*(v*3.0-1.0)*(s+t+u+v-1.0)*(-9.0/2.0);
  phi[13]= s*t*(s*3.0-1.0)*(9.0/2.0);
  phi[14]= s*t*(t*3.0-1.0)*(9.0/2.0);
  phi[15]= s*u*(s*3.0-1.0)*(9.0/2.0);
  phi[16]= s*u*(u*3.0-1.0)*(9.0/2.0);
  phi[17]= s*v*(s*3.0-1.0)*(9.0/2.0);
  phi[18]= s*v*(v*3.0-1.0)*(9.0/2.0);
  phi[19]= t*u*(t*3.0-1.0)*(9.0/2.0);
  phi[20]= t*u*(u*3.0-1.0)*(9.0/2.0);
  phi[21]= t*v*(t*3.0-1.0)*(9.0/2.0);
  phi[22]= t*v*(v*3.0-1.0)*(9.0/2.0);
  phi[23]= u*v*(u*3.0-1.0)*(9.0/2.0);
  phi[24]= u*v*(v*3.0-1.0)*(9.0/2.0);
  phi[25]= s*t*(s+t+u+v-1.0)*-2.7E1;
  phi[26]= s*u*(s+t+u+v-1.0)*-2.7E1;
  phi[27]= s*v*(s+t+u+v-1.0)*-2.7E1;
  phi[28]= t*u*(s+t+u+v-1.0)*-2.7E1;
  phi[29]= t*v*(s+t+u+v-1.0)*-2.7E1;
  phi[30]= u*v*(s+t+u+v-1.0)*-2.7E1;
  phi[31]= s*t*u*2.7E1;
  phi[32]= s*t*v*2.7E1;
  phi[33]= s*u*v*2.7E1;
  phi[34]= t*u*v*2.7E1;
}

template<>
void
Lagrange<Simplex,4,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  real_t* phis = phi;
  phis[0]= s*1.8E1+t*1.8E1+u*1.8E1+v*1.8E1-s*t*2.7E1-s*u*2.7E1-s*v*2.7E1-t*u*2.7E1-t*v*2.7E1-u*v*2.7E1-
             (s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-(u*u)*(2.7E1/2.0)-(v*v)*(2.7E1/2.0)-1.1E1/2.0;
  phis[1]= s*(-9.0/2.0)+(s*(s*1.8E1-9.0))/2.0+(s*s)*(9.0/2.0)+1.0;
  phis[2]= 0.0;
  phis[3]= 0.0;
  phis[4]= 0.0;
  phis[5]= s*(-4.5E1/2.0)-t*(4.5E1/2.0)-u*(4.5E1/2.0)-v*(4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+s*v*2.7E1+t*u*
             2.7E1+t*v*2.7E1+u*v*2.7E1+s*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0)+(s*s)*(2.7E1/2.0)+(t*t)*(2.7E1/2.0)+
             (u*u)*(2.7E1/2.0)+(v*v)*(2.7E1/2.0)+9.0;
  phis[6]= s*(s+t+u+v-1.0)*(-2.7E1/2.0)-s*(s*3.0-1.0)*(9.0/2.0)-(s*3.0-1.0)*(s+t+u+v-1.0)*(9.0/2.0);
  phis[7]= t*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phis[8]= t*(t*3.0-1.0)*(-9.0/2.0);
  phis[9]= u*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phis[10]= u*(u*3.0-1.0)*(-9.0/2.0);
  phis[11]= v*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phis[12]= v*(v*3.0-1.0)*(-9.0/2.0);
  phis[13]= t*(s*3.0-1.0)*(9.0/2.0)+s*t*(2.7E1/2.0);
  phis[14]= t*(t*3.0-1.0)*(9.0/2.0);
  phis[15]= u*(s*3.0-1.0)*(9.0/2.0)+s*u*(2.7E1/2.0);
  phis[16]= u*(u*3.0-1.0)*(9.0/2.0);
  phis[17]= v*(s*3.0-1.0)*(9.0/2.0)+s*v*(2.7E1/2.0);
  phis[18]= v*(v*3.0-1.0)*(9.0/2.0);
  phis[19]= 0.0;
  phis[20]= 0.0;
  phis[21]= 0.0;
  phis[22]= 0.0;
  phis[23]= 0.0;
  phis[24]= 0.0;
  phis[25]= t*(s+t+u+v-1.0)*-2.7E1-s*t*2.7E1;
  phis[26]= u*(s+t+u+v-1.0)*-2.7E1-s*u*2.7E1;
  phis[27]= v*(s+t+u+v-1.0)*-2.7E1-s*v*2.7E1;
  phis[28]= t*u*-2.7E1;
  phis[29]= t*v*-2.7E1;
  phis[30]= u*v*-2.7E1;
  phis[31]= t*u*2.7E1;
  phis[32]= t*v*2.7E1;
  phis[33]= u*v*2.7E1;
  phis[34]= 0.0;

  real_t* phit = phi + 35;
  phit[0]= s*1.8E1+t*1.8E1+u*1.8E1+v*1.8E1-s*t*2.7E1-s*u*2.7E1-s*v*2.7E1-t*u*2.7E1-t*v*2.7E1-u*v*2.7E1-
             (s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-(u*u)*(2.7E1/2.0)-(v*v)*(2.7E1/2.0)-1.1E1/2.0;
  phit[1]= 0.0;
  phit[2]= t*(-9.0/2.0)+(t*(t*1.8E1-9.0))/2.0+(t*t)*(9.0/2.0)+1.0;
  phit[3]= 0.0;
  phit[4]= 0.0;
  phit[5]= s*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phit[6]= s*(s*3.0-1.0)*(-9.0/2.0);
  phit[7]= s*(-4.5E1/2.0)-t*(4.5E1/2.0)-u*(4.5E1/2.0)-v*(4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+s*v*2.7E1+t*u*
             2.7E1+t*v*2.7E1+u*v*2.7E1+t*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0)+(s*s)*(2.7E1/2.0)+(t*t)*(2.7E1/2.0)+
             (u*u)*(2.7E1/2.0)+(v*v)*(2.7E1/2.0)+9.0;
  phit[8]= t*(s+t+u+v-1.0)*(-2.7E1/2.0)-t*(t*3.0-1.0)*(9.0/2.0)-(t*3.0-1.0)*(s+t+u+v-1.0)*(9.0/2.0);
  phit[9]= u*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phit[10]= u*(u*3.0-1.0)*(-9.0/2.0);
  phit[11]= v*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phit[12]= v*(v*3.0-1.0)*(-9.0/2.0);
  phit[13]= s*(s*3.0-1.0)*(9.0/2.0);
  phit[14]= s*(t*3.0-1.0)*(9.0/2.0)+s*t*(2.7E1/2.0);
  phit[15]= 0.0;
  phit[16]= 0.0;
  phit[17]= 0.0;
  phit[18]= 0.0;
  phit[19]= u*(t*3.0-1.0)*(9.0/2.0)+t*u*(2.7E1/2.0);
  phit[20]= u*(u*3.0-1.0)*(9.0/2.0);
  phit[21]= v*(t*3.0-1.0)*(9.0/2.0)+t*v*(2.7E1/2.0);
  phit[22]= v*(v*3.0-1.0)*(9.0/2.0);
  phit[23]= 0.0;
  phit[24]= 0.0;
  phit[25]= s*(s+t+u+v-1.0)*-2.7E1-s*t*2.7E1;
  phit[26]= s*u*-2.7E1;
  phit[27]= s*v*-2.7E1;
  phit[28]= u*(s+t+u+v-1.0)*-2.7E1-t*u*2.7E1;
  phit[29]= v*(s+t+u+v-1.0)*-2.7E1-t*v*2.7E1;
  phit[30]= u*v*-2.7E1;
  phit[31]= s*u*2.7E1;
  phit[32]= s*v*2.7E1;
  phit[33]= 0.0;
  phit[34]= u*v*2.7E1;

  real_t* phiu = phi + 70;
  phiu[0]= s*1.8E1+t*1.8E1+u*1.8E1+v*1.8E1-s*t*2.7E1-s*u*2.7E1-s*v*2.7E1-t*u*2.7E1-t*v*2.7E1-u*v*2.7E1-
             (s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-(u*u)*(2.7E1/2.0)-(v*v)*(2.7E1/2.0)-1.1E1/2.0;
  phiu[1]= 0.0;
  phiu[2]= 0.0;
  phiu[3]= u*(-9.0/2.0)+(u*(u*1.8E1-9.0))/2.0+(u*u)*(9.0/2.0)+1.0;
  phiu[4]= 0.0;
  phiu[5]= s*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiu[6]= s*(s*3.0-1.0)*(-9.0/2.0);
  phiu[7]= t*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiu[8]= t*(t*3.0-1.0)*(-9.0/2.0);
  phiu[9]= s*(-4.5E1/2.0)-t*(4.5E1/2.0)-u*(4.5E1/2.0)-v*(4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+s*v*2.7E1+t*u*
             2.7E1+t*v*2.7E1+u*v*2.7E1+u*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0)+(s*s)*(2.7E1/2.0)+(t*t)*(2.7E1/2.0)+
             (u*u)*(2.7E1/2.0)+(v*v)*(2.7E1/2.0)+9.0;
  phiu[10]= u*(s+t+u+v-1.0)*(-2.7E1/2.0)-u*(u*3.0-1.0)*(9.0/2.0)-(u*3.0-1.0)*(s+t+u+v-1.0)*(9.0/2.0);
  phiu[11]= v*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiu[12]= v*(v*3.0-1.0)*(-9.0/2.0);
  phiu[13]= 0.0;
  phiu[14]= 0.0;
  phiu[15]= s*(s*3.0-1.0)*(9.0/2.0);
  phiu[16]= s*(u*3.0-1.0)*(9.0/2.0)+s*u*(2.7E1/2.0);
  phiu[17]= 0.0;
  phiu[18]= 0.0;
  phiu[19]= t*(t*3.0-1.0)*(9.0/2.0);
  phiu[20]= t*(u*3.0-1.0)*(9.0/2.0)+t*u*(2.7E1/2.0);
  phiu[21]= 0.0;
  phiu[22]= 0.0;
  phiu[23]= v*(u*3.0-1.0)*(9.0/2.0)+u*v*(2.7E1/2.0);
  phiu[24]= v*(v*3.0-1.0)*(9.0/2.0);
  phiu[25]= s*t*-2.7E1;
  phiu[26]= s*(s+t+u+v-1.0)*-2.7E1-s*u*2.7E1;
  phiu[27]= s*v*-2.7E1;
  phiu[28]= t*(s+t+u+v-1.0)*-2.7E1-t*u*2.7E1;
  phiu[29]= t*v*-2.7E1;
  phiu[30]= v*(s+t+u+v-1.0)*-2.7E1-u*v*2.7E1;
  phiu[31]= s*t*2.7E1;
  phiu[32]= 0.0;
  phiu[33]= s*v*2.7E1;
  phiu[34]= t*v*2.7E1;

  real_t *phiv = phi + 105;
  phiv[0]= s*1.8E1+t*1.8E1+u*1.8E1+v*1.8E1-s*t*2.7E1-s*u*2.7E1-s*v*2.7E1-t*u*2.7E1-t*v*2.7E1-u*v*2.7E1-
             (s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-(u*u)*(2.7E1/2.0)-(v*v)*(2.7E1/2.0)-1.1E1/2.0;
  phiv[1]= 0.0;
  phiv[2]= 0.0;
  phiv[3]= 0.0;
  phiv[4]= v*(-9.0/2.0)+(v*(v*1.8E1-9.0))/2.0+(v*v)*(9.0/2.0)+1.0;
  phiv[5]= s*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiv[6]= s*(s*3.0-1.0)*(-9.0/2.0);
  phiv[7]= t*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiv[8]= t*(t*3.0-1.0)*(-9.0/2.0);
  phiv[9]= u*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0);
  phiv[10]= u*(u*3.0-1.0)*(-9.0/2.0);
  phiv[11]= s*(-4.5E1/2.0)-t*(4.5E1/2.0)-u*(4.5E1/2.0)-v*(4.5E1/2.0)+s*t*2.7E1+s*u*2.7E1+s*v*2.7E1+t*
             u*2.7E1+t*v*2.7E1+u*v*2.7E1+v*(s*6.0+t*6.0+u*6.0+v*6.0-5.0)*(9.0/2.0)+(s*s)*(2.7E1/2.0)+(t*t)*(2.7E1/
             2.0)+(u*u)*(2.7E1/2.0)+(v*v)*(2.7E1/2.0)+9.0;
  phiv[12]= v*(s+t+u+v-1.0)*(-2.7E1/2.0)-v*(v*3.0-1.0)*(9.0/2.0)-(v*3.0-1.0)*(s+t+u+v-1.0)*(9.0/2.0);
  phiv[13]= 0.0;
  phiv[14]= 0.0;
  phiv[15]= 0.0;
  phiv[16]= 0.0;
  phiv[17]= s*(s*3.0-1.0)*(9.0/2.0);
  phiv[18]= s*(v*3.0-1.0)*(9.0/2.0)+s*v*(2.7E1/2.0);
  phiv[19]= 0.0;
  phiv[20]= 0.0;
  phiv[21]= t*(t*3.0-1.0)*(9.0/2.0);
  phiv[22]= t*(v*3.0-1.0)*(9.0/2.0)+t*v*(2.7E1/2.0);
  phiv[23]= u*(u*3.0-1.0)*(9.0/2.0);
  phiv[24]= u*(v*3.0-1.0)*(9.0/2.0)+u*v*(2.7E1/2.0);
  phiv[25]= s*t*-2.7E1;
  phiv[26]= s*u*-2.7E1;
  phiv[27]= s*(s+t+u+v-1.0)*-2.7E1-s*v*2.7E1;
  phiv[28]= t*u*-2.7E1;
  phiv[29]= t*(s+t+u+v-1.0)*-2.7E1-t*v*2.7E1;
  phiv[30]= u*(s+t+u+v-1.0)*-2.7E1-u*v*2.7E1;
  phiv[31]= 0.0;
  phiv[32]= s*t*2.7E1;
  phiv[33]= s*u*2.7E1;
  phiv[34]= t*u*2.7E1;
}

template<>
void
Lagrange<Simplex,4,3>::hess( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  // phiss
  real_t *phiss = phi;
  phiss[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phiss[1]= s*2.7E1-9.0;
  phiss[2]= 0.0;
  phiss[3]= 0.0;
  phiss[4]= 0.0;
  phiss[5]= s*8.1E1+t*5.4E1+u*5.4E1+v*5.4E1-4.5E1;
  phiss[6]= s*-8.1E1-t*2.7E1-u*2.7E1-v*2.7E1+3.6E1;
  phiss[7]= t*2.7E1;
  phiss[8]= 0.0;
  phiss[9]= u*2.7E1;
  phiss[10]= 0.0;
  phiss[11]= v*2.7E1;
  phiss[12]= 0.0;
  phiss[13]= t*2.7E1;
  phiss[14]= 0.0;
  phiss[15]= u*2.7E1;
  phiss[16]= 0.0;
  phiss[17]= v*2.7E1;
  phiss[18]= 0.0;
  phiss[19]= 0.0;
  phiss[20]= 0.0;
  phiss[21]= 0.0;
  phiss[22]= 0.0;
  phiss[23]= 0.0;
  phiss[24]= 0.0;
  phiss[25]= t*-5.4E1;
  phiss[26]= u*-5.4E1;
  phiss[27]= v*-5.4E1;
  phiss[28]= 0.0;
  phiss[29]= 0.0;
  phiss[30]= 0.0;
  phiss[31]= 0.0;
  phiss[32]= 0.0;
  phiss[33]= 0.0;
  phiss[34]= 0.0;

  // phist
  real_t *phist = phi + 35;
  phist[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phist[1]= 0.0;
  phist[2]= 0.0;
  phist[3]= 0.0;
  phist[4]= 0.0;
  phist[5]= s*5.4E1+t*2.7E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phist[6]= s*-2.7E1+9.0/2.0;
  phist[7]= s*2.7E1+t*5.4E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phist[8]= t*-2.7E1+9.0/2.0;
  phist[9]= u*2.7E1;
  phist[10]= 0.0;
  phist[11]= v*2.7E1;
  phist[12]= 0.0;
  phist[13]= s*2.7E1-9.0/2.0;
  phist[14]= t*2.7E1-9.0/2.0;
  phist[15]= 0.0;
  phist[16]= 0.0;
  phist[17]= 0.0;
  phist[18]= 0.0;
  phist[19]= 0.0;
  phist[20]= 0.0;
  phist[21]= 0.0;
  phist[22]= 0.0;
  phist[23]= 0.0;
  phist[24]= 0.0;
  phist[25]= s*-5.4E1-t*5.4E1-u*2.7E1-v*2.7E1+2.7E1;
  phist[26]= u*-2.7E1;
  phist[27]= v*-2.7E1;
  phist[28]= u*-2.7E1;
  phist[29]= v*-2.7E1;
  phist[30]= 0.0;
  phist[31]= u*2.7E1;
  phist[32]= v*2.7E1;
  phist[33]= 0.0;
  phist[34]= 0.0;

  // phitt
  real_t* phitt = phi + 70;
  phitt[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phitt[1]= 0.0;
  phitt[2]= t*2.7E1-9.0;
  phitt[3]= 0.0;
  phitt[4]= 0.0;
  phitt[5]= s*2.7E1;
  phitt[6]= 0.0;
  phitt[7]= s*5.4E1+t*8.1E1+u*5.4E1+v*5.4E1-4.5E1;
  phitt[8]= s*-2.7E1-t*8.1E1-u*2.7E1-v*2.7E1+3.6E1;
  phitt[9]= u*2.7E1;
  phitt[10]= 0.0;
  phitt[11]= v*2.7E1;
  phitt[12]= 0.0;
  phitt[13]= 0.0;
  phitt[14]= s*2.7E1;
  phitt[15]= 0.0;
  phitt[16]= 0.0;
  phitt[17]= 0.0;
  phitt[18]= 0.0;
  phitt[19]= u*2.7E1;
  phitt[20]= 0.0;
  phitt[21]= v*2.7E1;
  phitt[22]= 0.0;
  phitt[23]= 0.0;
  phitt[24]= 0.0;
  phitt[25]= s*-5.4E1;
  phitt[26]= 0.0;
  phitt[27]= 0.0;
  phitt[28]= u*-5.4E1;
  phitt[29]= v*-5.4E1;
  phitt[30]= 0.0;
  phitt[31]= 0.0;
  phitt[32]= 0.0;
  phitt[33]= 0.0;
  phitt[34]= 0.0;

  // phisu
  real_t* phisu = phi + 105;
  phisu[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phisu[1]= 0.0;
  phisu[2]= 0.0;
  phisu[3]= 0.0;
  phisu[4]= 0.0;
  phisu[5]= s*5.4E1+t*2.7E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phisu[6]= s*-2.7E1+9.0/2.0;
  phisu[7]= t*2.7E1;
  phisu[8]= 0.0;
  phisu[9]= s*2.7E1+t*2.7E1+u*5.4E1+v*2.7E1-4.5E1/2.0;
  phisu[10]= u*-2.7E1+9.0/2.0;
  phisu[11]= v*2.7E1;
  phisu[12]= 0.0;
  phisu[13]= 0.0;
  phisu[14]= 0.0;
  phisu[15]= s*2.7E1-9.0/2.0;
  phisu[16]= u*2.7E1-9.0/2.0;
  phisu[17]= 0.0;
  phisu[18]= 0.0;
  phisu[19]= 0.0;
  phisu[20]= 0.0;
  phisu[21]= 0.0;
  phisu[22]= 0.0;
  phisu[23]= 0.0;
  phisu[24]= 0.0;
  phisu[25]= t*-2.7E1;
  phisu[26]= s*-5.4E1-t*2.7E1-u*5.4E1-v*2.7E1+2.7E1;
  phisu[27]= v*-2.7E1;
  phisu[28]= t*-2.7E1;
  phisu[29]= 0.0;
  phisu[30]= v*-2.7E1;
  phisu[31]= t*2.7E1;
  phisu[32]= 0.0;
  phisu[33]= v*2.7E1;
  phisu[34]= 0.0;

  // phitu
  real_t *phitu = phi + 140;
  phitu[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phitu[1]= 0.0;
  phitu[2]= 0.0;
  phitu[3]= 0.0;
  phitu[4]= 0.0;
  phitu[5]= s*2.7E1;
  phitu[6]= 0.0;
  phitu[7]= s*2.7E1+t*5.4E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phitu[8]= t*-2.7E1+9.0/2.0;
  phitu[9]= s*2.7E1+t*2.7E1+u*5.4E1+v*2.7E1-4.5E1/2.0;
  phitu[10]= u*-2.7E1+9.0/2.0;
  phitu[11]= v*2.7E1;
  phitu[12]= 0.0;
  phitu[13]= 0.0;
  phitu[14]= 0.0;
  phitu[15]= 0.0;
  phitu[16]= 0.0;
  phitu[17]= 0.0;
  phitu[18]= 0.0;
  phitu[19]= t*2.7E1-9.0/2.0;
  phitu[20]= u*2.7E1-9.0/2.0;
  phitu[21]= 0.0;
  phitu[22]= 0.0;
  phitu[23]= 0.0;
  phitu[24]= 0.0;
  phitu[25]= s*-2.7E1;
  phitu[26]= s*-2.7E1;
  phitu[27]= 0.0;
  phitu[28]= s*-2.7E1-t*5.4E1-u*5.4E1-v*2.7E1+2.7E1;
  phitu[29]= v*-2.7E1;
  phitu[30]= v*-2.7E1;
  phitu[31]= s*2.7E1;
  phitu[32]= 0.0;
  phitu[33]= 0.0;
  phitu[34]= v*2.7E1;

  // phiuu
  real_t *phiuu = phi + 175;
  phiuu[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phiuu[1]= 0.0;
  phiuu[2]= 0.0;
  phiuu[3]= u*2.7E1-9.0;
  phiuu[4]= 0.0;
  phiuu[5]= s*2.7E1;
  phiuu[6]= 0.0;
  phiuu[7]= t*2.7E1;
  phiuu[8]= 0.0;
  phiuu[9]= s*5.4E1+t*5.4E1+u*8.1E1+v*5.4E1-4.5E1;
  phiuu[10]= s*-2.7E1-t*2.7E1-u*8.1E1-v*2.7E1+3.6E1;
  phiuu[11]= v*2.7E1;
  phiuu[12]= 0.0;
  phiuu[13]= 0.0;
  phiuu[14]= 0.0;
  phiuu[15]= 0.0;
  phiuu[16]= s*2.7E1;
  phiuu[17]= 0.0;
  phiuu[18]= 0.0;
  phiuu[19]= 0.0;
  phiuu[20]= t*2.7E1;
  phiuu[21]= 0.0;
  phiuu[22]= 0.0;
  phiuu[23]= v*2.7E1;
  phiuu[24]= 0.0;
  phiuu[25]= 0.0;
  phiuu[26]= s*-5.4E1;
  phiuu[27]= 0.0;
  phiuu[28]= t*-5.4E1;
  phiuu[29]= 0.0;
  phiuu[30]= v*-5.4E1;
  phiuu[31]= 0.0;
  phiuu[32]= 0.0;
  phiuu[33]= 0.0;
  phiuu[34]= 0.0;

  // phisv
  real_t *phisv = phi + 210;
  phisv[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phisv[1]= 0.0;
  phisv[2]= 0.0;
  phisv[3]= 0.0;
  phisv[4]= 0.0;
  phisv[5]= s*5.4E1+t*2.7E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phisv[6]= s*-2.7E1+9.0/2.0;
  phisv[7]= t*2.7E1;
  phisv[8]= 0.0;
  phisv[9]= u*2.7E1;
  phisv[10]= 0.0;
  phisv[11]= s*2.7E1+t*2.7E1+u*2.7E1+v*5.4E1-4.5E1/2.0;
  phisv[12]= v*-2.7E1+9.0/2.0;
  phisv[13]= 0.0;
  phisv[14]= 0.0;
  phisv[15]= 0.0;
  phisv[16]= 0.0;
  phisv[17]= s*2.7E1-9.0/2.0;
  phisv[18]= v*2.7E1-9.0/2.0;
  phisv[19]= 0.0;
  phisv[20]= 0.0;
  phisv[21]= 0.0;
  phisv[22]= 0.0;
  phisv[23]= 0.0;
  phisv[24]= 0.0;
  phisv[25]= t*-2.7E1;
  phisv[26]= u*-2.7E1;
  phisv[27]= s*-5.4E1-t*2.7E1-u*2.7E1-v*5.4E1+2.7E1;
  phisv[28]= 0.0;
  phisv[29]= t*-2.7E1;
  phisv[30]= u*-2.7E1;
  phisv[31]= 0.0;
  phisv[32]= t*2.7E1;
  phisv[33]= u*2.7E1;
  phisv[34]= 0.0;

  // phitv
  real_t *phitv = phi + 245;
  phitv[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phitv[1]= 0.0;
  phitv[2]= 0.0;
  phitv[3]= 0.0;
  phitv[4]= 0.0;
  phitv[5]= s*2.7E1;
  phitv[6]= 0.0;
  phitv[7]= s*2.7E1+t*5.4E1+u*2.7E1+v*2.7E1-4.5E1/2.0;
  phitv[8]= t*-2.7E1+9.0/2.0;
  phitv[9]= u*2.7E1;
  phitv[10]= 0.0;
  phitv[11]= s*2.7E1+t*2.7E1+u*2.7E1+v*5.4E1-4.5E1/2.0;
  phitv[12]= v*-2.7E1+9.0/2.0;
  phitv[13]= 0.0;
  phitv[14]= 0.0;
  phitv[15]= 0.0;
  phitv[16]= 0.0;
  phitv[17]= 0.0;
  phitv[18]= 0.0;
  phitv[19]= 0.0;
  phitv[20]= 0.0;
  phitv[21]= t*2.7E1-9.0/2.0;
  phitv[22]= v*2.7E1-9.0/2.0;
  phitv[23]= 0.0;
  phitv[24]= 0.0;
  phitv[25]= s*-2.7E1;
  phitv[26]= 0.0;
  phitv[27]= s*-2.7E1;
  phitv[28]= u*-2.7E1;
  phitv[29]= s*-2.7E1-t*5.4E1-u*2.7E1-v*5.4E1+2.7E1;
  phitv[30]= u*-2.7E1;
  phitv[31]= 0.0;
  phitv[32]= s*2.7E1;
  phitv[33]= 0.0;
  phitv[34]= u*2.7E1;

  // phiuv
  real_t *phiuv = phi + 280;
  phiuv[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phiuv[1]= 0.0;
  phiuv[2]= 0.0;
  phiuv[3]= 0.0;
  phiuv[4]= 0.0;
  phiuv[5]= s*2.7E1;
  phiuv[6]= 0.0;
  phiuv[7]= t*2.7E1;
  phiuv[8]= 0.0;
  phiuv[9]= s*2.7E1+t*2.7E1+u*5.4E1+v*2.7E1-4.5E1/2.0;
  phiuv[10]= u*-2.7E1+9.0/2.0;
  phiuv[11]= s*2.7E1+t*2.7E1+u*2.7E1+v*5.4E1-4.5E1/2.0;
  phiuv[12]= v*-2.7E1+9.0/2.0;
  phiuv[13]= 0.0;
  phiuv[14]= 0.0;
  phiuv[15]= 0.0;
  phiuv[16]= 0.0;
  phiuv[17]= 0.0;
  phiuv[18]= 0.0;
  phiuv[19]= 0.0;
  phiuv[20]= 0.0;
  phiuv[21]= 0.0;
  phiuv[22]= 0.0;
  phiuv[23]= u*2.7E1-9.0/2.0;
  phiuv[24]= v*2.7E1-9.0/2.0;
  phiuv[25]= 0.0;
  phiuv[26]= s*-2.7E1;
  phiuv[27]= s*-2.7E1;
  phiuv[28]= t*-2.7E1;
  phiuv[29]= t*-2.7E1;
  phiuv[30]= s*-2.7E1-t*2.7E1-u*5.4E1-v*5.4E1+2.7E1;
  phiuv[31]= 0.0;
  phiuv[32]= 0.0;
  phiuv[33]= s*2.7E1;
  phiuv[34]= t*2.7E1;

  // phivv
  real_t *phivv = phi + 315;
  phivv[0]= s*-2.7E1-t*2.7E1-u*2.7E1-v*2.7E1+1.8E1;
  phivv[1]= 0.0;
  phivv[2]= 0.0;
  phivv[3]= 0.0;
  phivv[4]= v*2.7E1-9.0;
  phivv[5]= s*2.7E1;
  phivv[6]= 0.0;
  phivv[7]= t*2.7E1;
  phivv[8]= 0.0;
  phivv[9]= u*2.7E1;
  phivv[10]= 0.0;
  phivv[11]= s*5.4E1+t*5.4E1+u*5.4E1+v*8.1E1-4.5E1;
  phivv[12]= s*-2.7E1-t*2.7E1-u*2.7E1-v*8.1E1+3.6E1;
  phivv[13]= 0.0;
  phivv[14]= 0.0;
  phivv[15]= 0.0;
  phivv[16]= 0.0;
  phivv[17]= 0.0;
  phivv[18]= s*2.7E1;
  phivv[19]= 0.0;
  phivv[20]= 0.0;
  phivv[21]= 0.0;
  phivv[22]= t*2.7E1;
  phivv[23]= 0.0;
  phivv[24]= u*2.7E1;
  phivv[25]= 0.0;
  phivv[26]= 0.0;
  phivv[27]= s*-5.4E1;
  phivv[28]= 0.0;
  phivv[29]= t*-5.4E1;
  phivv[30]= u*-5.4E1;
  phivv[31]= 0.0;
  phivv[32]= 0.0;
  phivv[33]= 0.0;
  phivv[34]= 0.0;
}

/*
 * pentatopes, p = 4
*/
template<>
void
Lagrange<Simplex,4,4>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,4,4>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,4,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * pentatopes, p = 5
*/
template<>
void
Lagrange<Simplex,4,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,4,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Lagrange<Simplex,4,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
