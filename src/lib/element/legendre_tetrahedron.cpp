#include "common/error.h"
#include "element/basis.h"

#include <cmath>

namespace avro
{

/*
 * tetrahedron, p = 0
*/

template<>
void
Legendre<Simplex,3,0>::eval( const real_t* x , real_t* phi ) {
  phi[0] = 1;
}

template<>
void
Legendre<Simplex,3,0>::grad( const real_t* x , real_t* phi ) {
  phi[0] = 0;
  phi[0] = 0;
  phi[0] = 0;
}

template<>
void
Legendre<Simplex,3,0>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 6; i++)
    phi[i] = 0;
}

/*
 * tetrahedron, p = 1
*/

template<>
void
Legendre<Simplex,3,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;

  phi[0] = 1;
  phi[1] = (sqrt(5/3.))*(-1+(4)*(s));
  phi[2] = (sqrt(10/3.))*(-1+s+(3)*(t));
  phi[3] = (sqrt(10))*(-v+(2)*(u));
}

template<>
void
Legendre<Simplex,3,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (4)*(sqrt(5/3.));
  phis[2] = sqrt(10/3.);
  phis[3] = sqrt(10);

  real_t* phit = phi + 4;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = sqrt(30);
  phit[3] = sqrt(10);

  real_t* phiu = phi + 8;
  phiu[0] = 0;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = (2)*(sqrt(10));
}

template<>
void
Legendre<Simplex,3,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 24; i++)
    phi[i] = 0;
}

/*
 * pentatopes, p = 2
*/
template<>
void
Legendre<Simplex,3,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;
  real_t v2 = v*v;

  real_t t2 = t*t;

  real_t u2 = u*u;

  phi[0] = 1;
  phi[1] = (sqrt(5/3.))*(-1+(4)*(s));
  phi[2] = (sqrt(10/3.))*(-1+s+(3)*(t));
  phi[3] = (sqrt(10))*(-v+(2)*(u));
  phi[4] = (sqrt(7/3.))*(1+(5)*((s)*(-2+(3)*(s))));
  phi[5] = (sqrt(14/3.))*((-1+(6)*(s))*(-1+s+(3)*(t)));
  phi[6] = (sqrt(14))*((-1+(6)*(s))*(-v+(2)*(u)));
  phi[7] = (sqrt(7))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*(t2));
  phi[8] = (sqrt(21))*((-1+s+(5)*(t))*(-v+(2)*(u)));
  phi[9] = (sqrt(35))*(v2+(6)*((-v)*(u))+(6)*(u2));
}

template<>
void
Legendre<Simplex,3,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (4)*(sqrt(5/3.));
  phis[2] = sqrt(10/3.);
  phis[3] = sqrt(10);
  phis[4] = (10)*((sqrt(7/3.))*(-1+(3)*(s)));
  phis[5] = (sqrt(14/3.))*(-7+(12)*(s)+(18)*(t));
  phis[6] = (sqrt(14))*(-7+(12)*(s)+(6)*(t)+(12)*(u));
  phis[7] = (2)*((sqrt(7))*(-1+s+(4)*(t)));
  phis[8] = (2)*((sqrt(21))*(-1+s+(3)*(t)+u));
  phis[9] = (2)*((sqrt(35))*(-v+(3)*(u)));

  real_t* phit = phi + 10;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = sqrt(30);
  phit[3] = sqrt(10);
  phit[4] = 0;
  phit[5] = (sqrt(42))*(-1+(6)*(s));
  phit[6] = (sqrt(14))*(-1+(6)*(s));
  phit[7] = (4)*((sqrt(7))*(-2+(2)*(s)+(5)*(t)));
  phit[8] = (2)*((sqrt(21))*(-3+(3)*(s)+(5)*(t)+(5)*(u)));
  phit[9] = (2)*((sqrt(35))*(-v+(3)*(u)));

  real_t* phiu = phi + 20;
  phiu[0] = 0;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = (2)*(sqrt(10));
  phiu[4] = 0;
  phiu[5] = 0;
  phiu[6] = (2)*((sqrt(14))*(-1+(6)*(s)));
  phiu[7] = 0;
  phiu[8] = (2)*((sqrt(21))*(-1+s+(5)*(t)));
  phiu[9] = (6)*((sqrt(35))*(-v+(2)*(u)));

}

template<>
void
Legendre<Simplex,3,2>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * tetrahedron, p = 3
*/
template<>
void
Legendre<Simplex,3,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;
  real_t v2 = v*v;
  real_t v3 = v2*v;

  real_t t2 = t*t;
  real_t t3 = t2*t;

  real_t u2 = u*u;
  real_t u3 = u2*u;

  phi[0] = 1;
  phi[1] = (sqrt(5/3.))*(-1+(4)*(s));
  phi[2] = (sqrt(10/3.))*(-1+s+(3)*(t));
  phi[3] = (sqrt(10))*(-v+(2)*(u));
  phi[4] = (sqrt(7/3.))*(1+(5)*((s)*(-2+(3)*(s))));
  phi[5] = (sqrt(14/3.))*((-1+(6)*(s))*(-1+s+(3)*(t)));
  phi[6] = (sqrt(14))*((-1+(6)*(s))*(-v+(2)*(u)));
  phi[7] = (sqrt(7))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*(t2));
  phi[8] = (sqrt(21))*((-1+s+(5)*(t))*(-v+(2)*(u)));
  phi[9] = (sqrt(35))*(v2+(6)*((-v)*(u))+(6)*(u2));
  phi[10] = (sqrt(3))*(-1+(s)*(18+(7)*((s)*(-9+(8)*(s)))));
  phi[11] = (sqrt(6))*((1+(14)*((s)*(-1+(2)*(s))))*(-1+s+(3)*(t)));
  phi[12] = (3)*((sqrt(2))*((1+(14)*((s)*(-1+(2)*(s))))*(-v+(2)*(u))));
  phi[13] = (3)*((-1+(8)*(s))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*(t2)));
  phi[14] = (3)*((sqrt(3))*((-1+(8)*(s))*((-1+s+(5)*(t))*(-v+(2)*(u)))));
  phi[15] = (3)*((sqrt(5))*((-1+(8)*(s))*(v2+(6)*((-v)*(u))+(6)*(u2))));
  phi[16] = (2)*((sqrt(3))*(pow(-1+s,3)+(15)*((pow(-1+s,2))*(t))+(45)*((-1+s)*(t2))+(35)*(t3)));
  phi[17] = (6)*((pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*(t2))*(-v+(2)*(u)));
  phi[18] = (2)*((sqrt(15))*((-1+s+(7)*(t))*(v2+(6)*((-v)*(u))+(6)*(u2))));
  phi[19] = (2)*((sqrt(21))*((-v3)+(12)*(v2*(u))+(30)*((-v)*(u2))+(20)*(u3)));
}

template<>
void
Legendre<Simplex,3,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;
  real_t v2 = v*v;

  real_t s2 = s*s;
  real_t t2 = t*t;
  real_t u2 = u*u;

  real_t* phis = phi;

  phis[0] = 0;
  phis[1] = (4)*(sqrt(5/3.));
  phis[2] = sqrt(10/3.);
  phis[3] = sqrt(10);
  phis[4] = (10)*((sqrt(7/3.))*(-1+(3)*(s)));
  phis[5] = (sqrt(14/3.))*(-7+(12)*(s)+(18)*(t));
  phis[6] = (sqrt(14))*(-7+(12)*(s)+(6)*(t)+(12)*(u));
  phis[7] = (2)*((sqrt(7))*(-1+s+(4)*(t)));
  phis[8] = (2)*((sqrt(21))*(-1+s+(3)*(t)+u));
  phis[9] = (2)*((sqrt(35))*(-v+(3)*(u)));
  phis[10] = (6)*((sqrt(3))*(3+(7)*((s)*(-3+(4)*(s)))));
  phis[11] = (3)*((sqrt(6))*(5+(-14)*(t)+(28)*((s)*(-1+s+(2)*(t)))));
  phis[12] = (3)*((sqrt(2))*(15+(-14)*(t)+(-28)*(u)+(28)*((s)*(-3+(3)*(s)+(2)*(t)+(4)*(u)))));
  phis[13] = (6)*(5+(12)*(s2)+(4)*((t)*(-9+(10)*(t)))+(s)*(-17+(64)*(t)));
  phis[14] = (6)*((sqrt(3))*(5+(12)*(s2)+(-9)*(u)+(s)*(-17+(48)*(t)+(16)*(u))+(t)*(-27+(20)*(t)+(40)*(u))));
  phis[15] = (6)*((sqrt(5))*(5+(12)*(s2)+(4)*(t2)+(3)*((u)*(-9+(8)*(u)))+(3)*((t)*(-3+(8)*(u)))+(s)*(-17+(16)*(t)+(48)*(u))));
  phis[16] = (6)*((sqrt(3))*(pow(-1+s,2)+(10)*((-1+s)*(t))+(15)*(t2)));
  phis[17] = (6)*(3+(3)*(s2)+(-4)*(u)+(s)*(-6+(26)*(t)+(4)*(u))+(t)*(-26+(33)*(t)+(24)*(u)));
  phis[18] = (6)*((sqrt(15))*(1+s2+(5)*(t2)+(2)*((-2+u)*(u))+(s)*(-2+(6)*(t)+(4)*(u))+(2)*((t)*(-3+(8)*(u)))));
  phis[19] = (6)*((sqrt(21))*(v2+(8)*((-v)*(u))+(10)*(u2)));

  real_t* phit = phi + 20;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = sqrt(30);
  phit[3] = sqrt(10);
  phit[4] = 0;
  phit[5] = (sqrt(42))*(-1+(6)*(s));
  phit[6] = (sqrt(14))*(-1+(6)*(s));
  phit[7] = (4)*((sqrt(7))*(-2+(2)*(s)+(5)*(t)));
  phit[8] = (2)*((sqrt(21))*(-3+(3)*(s)+(5)*(t)+(5)*(u)));
  phit[9] = (2)*((sqrt(35))*(-v+(3)*(u)));
  phit[10] = 0;
  phit[11] = (3)*((sqrt(6))*(1+(14)*((s)*(-1+(2)*(s)))));
  phit[12] = (3)*((sqrt(2))*(1+(14)*((s)*(-1+(2)*(s)))));
  phit[13] = (12)*((-1+(8)*(s))*(-2+(2)*(s)+(5)*(t)));
  phit[14] = (6)*((sqrt(3))*((-1+(8)*(s))*(-3+(3)*(s)+(5)*(t)+(5)*(u))));
  phit[15] = (6)*((sqrt(5))*((-1+(8)*(s))*(-v+(3)*(u))));
  phit[16] = (30)*((sqrt(3))*(pow(-1+s,2)+(6)*((-1+s)*(t))+(7)*(t2)));
  phit[17] = (6)*(13+(13)*(s2)+(-24)*(u)+(s)*(-26+(66)*(t)+(24)*(u))+(3)*((t)*(-22+(21)*(t)+(28)*(u))));
  phit[18] = (6)*((sqrt(15))*(3+(3)*(s2)+(7)*(t2)+(2)*((u)*(-8+(7)*(u)))+(2)*((s)*(-3+(5)*(t)+(8)*(u)))+(2)*((t)*(-5+(14)*(u)))));
  phit[19] = (6)*((sqrt(21))*(v2+(8)*((-v)*(u))+(10)*(u2)));

  real_t* phiu = phi + 40;
  phiu[0] = 0;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = (2)*(sqrt(10));
  phiu[4] = 0;
  phiu[5] = 0;
  phiu[6] = (2)*((sqrt(14))*(-1+(6)*(s)));
  phiu[7] = 0;
  phiu[8] = (2)*((sqrt(21))*(-1+s+(5)*(t)));
  phiu[9] = (6)*((sqrt(35))*(-v+(2)*(u)));
  phiu[10] = 0;
  phiu[11] = 0;
  phiu[12] = (6)*((sqrt(2))*(1+(14)*((s)*(-1+(2)*(s)))));
  phiu[13] = 0;
  phiu[14] = (6)*((sqrt(3))*((-1+(8)*(s))*(-1+s+(5)*(t))));
  phiu[15] = (18)*((sqrt(5))*((-1+(8)*(s))*(-v+(2)*(u))));
  phiu[16] = 0;
  phiu[17] = (12)*(pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*(t2));
  phiu[18] = (12)*((sqrt(15))*((-1+s+(7)*(t))*(-v+(2)*(u))));
  phiu[19] = (24)*((sqrt(21))*(v2+(5)*((-v)*(u))+(5)*(u2)));
}

template<>
void
Legendre<Simplex,3,3>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * tetrahedra, p = 4
*/
template<>
void
Legendre<Simplex,3,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;
  real_t v2 = v*v;
  real_t v3 = v2*v;

  real_t t2 = t*t;
  real_t t3 = t2*t;

  real_t u2 = u*u;
  real_t u3 = u2*u;

  phi[0] = 1;
  phi[1] = (sqrt(5/3.))*(-1+(4)*(s));
  phi[2] = (sqrt(10/3.))*(-1+s+(3)*(t));
  phi[3] = (sqrt(10))*(-v+(2)*(u));
  phi[4] = (sqrt(7/3.))*(1+(5)*((s)*(-2+(3)*(s))));
  phi[5] = (sqrt(14/3.))*((-1+(6)*(s))*(-1+s+(3)*(t)));
  phi[6] = (sqrt(14))*((-1+(6)*(s))*(-v+(2)*(u)));
  phi[7] = (sqrt(7))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*(t2));
  phi[8] = (sqrt(21))*((-1+s+(5)*(t))*(-v+(2)*(u)));
  phi[9] = (sqrt(35))*(v2+(6)*((-v)*(u))+(6)*(u2));
  phi[10] = (sqrt(3))*(-1+(s)*(18+(7)*((s)*(-9+(8)*(s)))));
  phi[11] = (sqrt(6))*((1+(14)*((s)*(-1+(2)*(s))))*(-1+s+(3)*(t)));
  phi[12] = (3)*((sqrt(2))*((1+(14)*((s)*(-1+(2)*(s))))*(-v+(2)*(u))));
  phi[13] = (3)*((-1+(8)*(s))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*(t2)));
  phi[14] = (3)*((sqrt(3))*((-1+(8)*(s))*((-1+s+(5)*(t))*(-v+(2)*(u)))));
  phi[15] = (3)*((sqrt(5))*((-1+(8)*(s))*(v2+(6)*((-v)*(u))+(6)*(u2))));
  phi[16] = (2)*((sqrt(3))*(pow(-1+s,3)+(15)*((pow(-1+s,2))*(t))+(45)*((-1+s)*(t2))+(35)*(t3)));
  phi[17] = (6)*((pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*(t2))*(-v+(2)*(u)));
  phi[18] = (2)*((sqrt(15))*((-1+s+(7)*(t))*(v2+(6)*((-v)*(u))+(6)*(u2))));
  phi[19] = (2)*((sqrt(21))*((-v3)+(12)*(v2*(u))+(30)*((-v)*(u2))+(20)*(u3)));
  phi[20] = (sqrt(11/3.))*(1+(14)*((s)*(-2+(3)*((s)*(4+(s)*(-8+(5)*(s)))))));
  phi[21] = (sqrt(22/3.))*((-1+(12)*((s)*(2+(s)*(-9+(10)*(s)))))*(-1+s+(3)*(t)));
  phi[22] = (sqrt(22))*((-1+(12)*((s)*(2+(s)*(-9+(10)*(s)))))*(-v+(2)*(u)));
  phi[23] = (sqrt(11))*((1+(9)*((s)*(-2+(5)*(s))))*(pow(-1+s,2)+(8)*((-1+s)*(t))+(10)*t2));
  phi[24] = (sqrt(33))*((1+(9)*((s)*(-2+(5)*(s))))*((-1+s+(5)*(t))*(-v+(2)*(u))));
  phi[25] = (sqrt(55))*((1+(9)*((s)*(-2+(5)*(s))))*(v2+(6)*((-v)*(u))+(6)*u2));
  phi[26] = (2)*((sqrt(11/3.))*((-1+(10)*(s))*(pow(-1+s,3)+(15)*((pow(-1+s,2))*(t))+(45)*((-1+s)*t2)+(35)*t3)));
  phi[27] = (2)*((sqrt(11))*((-1+(10)*(s))*((pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*t2)*(-v+(2)*(u)))));
  phi[28] = (2)*((sqrt(55/3.))*((-1+(10)*(s))*((-1+s+(7)*(t))*(v2+(6)*((-v)*(u))+(6)*u2))));
  phi[29] = (2)*((sqrt(77/3.))*((-1+(10)*(s))*((-v3)+(12)*(v2*(u))+(30)*((-v)*u2)+(20)*u3)));
  phi[30] = (sqrt(55/3.))*(pow(-1+s,4)+(24)*((pow(-1+s,3))*(t))+(126)*((pow(-1+s,2))*t2)+(224)*((-1+s)*t3)+(126)*(pow(t,4)));
  phi[31] = (sqrt(55))*((pow(-1+s,3)+(21)*((pow(-1+s,2))*(t))+(84)*((-1+s)*t2)+(84)*t3)*(-v+(2)*(u)));
  phi[32] = (5)*((sqrt(11/3.))*((pow(-1+s,2)+(16)*((-1+s)*(t))+(36)*t2)*(v2+(6)*((-v)*(u))+(6)*u2)));
  phi[33] = (sqrt(385/3.))*((-1+s+(9)*(t))*((-v3)+(12)*(v2*(u))+(30)*((-v)*u2)+(20)*u3));
  phi[34] = (sqrt(165))*(pow(-v,4)+(20)*((-v3)*(u))+(90)*(v2*u2)+(140)*((-v)*u3)+(70)*(pow(u,4)));
}

template<>
void
Legendre<Simplex,3,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];

  real_t v = 1 - s - t;

  real_t s2 = s*s;
  real_t t2 = t*t;
  real_t u2 = u*u;
  real_t v2 = v*v;

  real_t s3 = s2*s;
  real_t t3 = t2*t;
  real_t u3 = u2*u;
  real_t v3 = v2*v;

  // phis
  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (4)*(sqrt(5/3.));
  phis[2] = sqrt(10/3.);
  phis[3] = sqrt(10);
  phis[4] = (10)*((sqrt(7/3.))*(-1+(3)*(s)));
  phis[5] = (sqrt(14/3.))*(-7+(12)*(s)+(18)*(t));
  phis[6] = (sqrt(14))*(-7+(12)*(s)+(6)*(t)+(12)*(u));
  phis[7] = (2)*((sqrt(7))*(-1+s+(4)*(t)));
  phis[8] = (2)*((sqrt(21))*(-1+s+(3)*(t)+u));
  phis[9] = (2)*((sqrt(35))*(-v+(3)*(u)));
  phis[10] = (6)*((sqrt(3))*(3+(7)*((s)*(-3+(4)*(s)))));
  phis[11] = (3)*((sqrt(6))*(5+(-14)*(t)+(28)*((s)*(-1+s+(2)*(t)))));
  phis[12] = (3)*((sqrt(2))*(15+(-14)*(t)+(-28)*(u)+(28)*((s)*(-3+(3)*(s)+(2)*(t)+(4)*(u)))));
  phis[13] = (6)*(5+(12)*(s2)+(4)*((t)*(-9+(10)*(t)))+(s)*(-17+(64)*(t)));
  phis[14] = (6)*((sqrt(3))*(5+(12)*(s2)+(-9)*(u)+(s)*(-17+(48)*(t)+(16)*(u))+(t)*(-27+(20)*(t)+(40)*(u))));
  phis[15] = (6)*((sqrt(5))*(5+(12)*(s2)+(4)*(t2)+(3)*((u)*(-9+(8)*(u)))+(3)*((t)*(-3+(8)*(u)))+(s)*(-17+(16)*(t)+(48)*(u))));
  phis[16] = (6)*((sqrt(3))*(pow(-1+s,2)+(10)*((-1+s)*(t))+(15)*(t2)));
  phis[17] = (6)*(3+(3)*(s2)+(-4)*(u)+(s)*(-6+(26)*(t)+(4)*(u))+(t)*(-26+(33)*(t)+(24)*(u)));
  phis[18] = (6)*((sqrt(15))*(1+s2+(5)*(t2)+(2)*((-2+u)*(u))+(s)*(-2+(6)*(t)+(4)*(u))+(2)*((t)*(-3+(8)*(u)))));
  phis[19] = (6)*((sqrt(21))*(v2+(8)*((-v)*(u))+(10)*(u2)));
  phis[20] = (28)*((sqrt(11/3.))*(-1+(6)*((s)*(2+(s)*(-6+(5)*(s))))));
  phis[21] = (sqrt(22/3.))*(-25+(72)*(t)+(12)*((s)*(22+(-54)*(t)+(s)*(-57+(40)*(s)+(90)*(t)))));
  phis[22] = (sqrt(22))*(-1+(24)*(s)+(-108)*s2+(120)*s3+(24)*((1+(3)*((s)*(-3+(5)*(s))))*(-v+(2)*(u))));
  phis[23] = (4)*((sqrt(11))*(-5+(s)*(41+(9)*((s)*(-9+(5)*(s))))+(38)*(t)+(18)*((s)*((-14+(15)*(s))*(t)))+(45)*((-1+(5)*(s))*t2)));
  phis[24] = (sqrt(33))*((1+(9)*((s)*(-2+(5)*(s))))*(-1+s+(5)*(t))+(1+(9)*((s)*(-2+(5)*(s))))*(-v+(2)*(u))
             +(18)*((-1+(5)*(s))*((-1+s+(5)*(t))*(-v+(2)*(u)))));
  phis[25] = (sqrt(55))*((2)*((1+(9)*((s)*(-2+(5)*(s))))*(-v+(3)*(u)))+(-18+(90)*(s))*(v2+(6)*((-v)*(u))+(6)*u2));
  phis[26] = (2)*((sqrt(11/3.))*((pow(-1+s,2))*(-13+(40)*(s))+(90)*((-1+s)*((-2+(5)*(s))*(t)))+(45)*((-11+(20)*(s))*t2)+(350)*t3));
  phis[27] = (2)*((sqrt(11))*(-13+(40)*s3+(24)*(u)+s2*(-93+(390)*(t)+(60)*(u))+(6)*((s)*(11+(-14)*(u)+(t)*(-91+(110)*(t)+(80)*(u))))
             +(3)*((t)*(52+(-88)*(u)+(t)*(-121+(70)*(t)+(140)*(u))))));
  phis[28] = (2)*((sqrt(55/3.))*(-13+(40)*s3+(70)*t3+(72)*(u)+(-66)*u2+(15)*(t2*(-11+(28)*(u)))+(12)*((t)*((-1+u)*(-9+(35)*(u))))
             +(3)*(s2*(-31+(90)*(t)+(60)*(u)))+(6)*((s)*(11+(t)*(-63+(50)*(t))+(-42)*(u)+(160)*((t)*(u))+(20)*u2))));
  phis[29] = (2)*((sqrt(77/3.))*(v2*(-13+(40)*(s)+(10)*(t))+(24)*((-v)*((-6+(15)*(s)+(5)*(t))*(u)))+(30)*((-11+(20)*(s)+(10)*(t))*u2)+(200)*u3));
  phis[30] = (4)*((sqrt(55/3.))*(pow(-1+s,3)+(18)*((pow(-1+s,2))*(t))+(63)*((-1+s)*t2)+(56)*t3));
  phis[31] = (2)*((sqrt(55))*(-2+(2)*s3+(3)*(u)+(3)*(s2*(-2+(11)*(t)+u))+(3)*((t)*(11+(-14)*(u)+(7)*((t)*(-5+(4)*(t)+(4)*(u)))))
             +(3)*((s)*(2+(-2)*(u)+(t)*(-22+(35)*(t)+(14)*(u))))));
  phis[32] = (5)*((sqrt(11/3.))*((2)*((pow(-1+s,2)+(16)*((-1+s)*(t))+(36)*t2)*(-v+(3)*(u)))+(2)*((-1+s+(8)*(t))*(v2+(6)*((-v)*(u))+(6)*u2))));
  phis[33] = (4)*((sqrt(385/3.))*(v2*(-1+s+(7)*(t))+(3)*((-v)*((-3+(3)*(s)+(19)*(t))*(u)))+(15)*((-1+s+(5)*(t))*u2)+(5)*u3));
  phis[34] = (4)*((sqrt(165))*((-v3)+(15)*(v2*(u))+(45)*((-v)*u2)+(35)*u3));

  // phit
  real_t *phit = phi + 35;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = sqrt(30);
  phit[3] = sqrt(10);
  phit[4] = 0;
  phit[5] = (sqrt(42))*(-1+(6)*(s));
  phit[6] = (sqrt(14))*(-1+(6)*(s));
  phit[7] = (4)*((sqrt(7))*(-2+(2)*(s)+(5)*(t)));
  phit[8] = (2)*((sqrt(21))*(-3+(3)*(s)+(5)*(t)+(5)*(u)));
  phit[9] = (2)*((sqrt(35))*(-v+(3)*(u)));
  phit[10] = 0;
  phit[11] = (3)*((sqrt(6))*(1+(14)*((s)*(-1+(2)*(s)))));
  phit[12] = (3)*((sqrt(2))*(1+(14)*((s)*(-1+(2)*(s)))));
  phit[13] = (12)*((-1+(8)*(s))*(-2+(2)*(s)+(5)*(t)));
  phit[14] = (6)*((sqrt(3))*((-1+(8)*(s))*(-3+(3)*(s)+(5)*(t)+(5)*(u))));
  phit[15] = (6)*((sqrt(5))*((-1+(8)*(s))*(-v+(3)*(u))));
  phit[16] = (30)*((sqrt(3))*(pow(-1+s,2)+(6)*((-1+s)*(t))+(7)*(t2)));
  phit[17] = (6)*(13+(13)*(s2)+(-24)*(u)+(s)*(-26+(66)*(t)+(24)*(u))+(3)*((t)*(-22+(21)*(t)+(28)*(u))));
  phit[18] = (6)*((sqrt(15))*(3+(3)*(s2)+(7)*(t2)+(2)*((u)*(-8+(7)*(u)))+(2)*((s)*(-3+(5)*(t)+(8)*(u)))+(2)*((t)*(-5+(14)*(u)))));
  phit[19] = (6)*((sqrt(21))*(v2+(8)*((-v)*(u))+(10)*(u2)));
  phit[20] = 0;
  phit[21] = (sqrt(66))*(-1+(12)*((s)*(2+(s)*(-9+(10)*(s)))));
  phit[22] = (sqrt(22))*(-1+(12)*((s)*(2+(s)*(-9+(10)*(s)))));
  phit[23] = (4)*((sqrt(11))*((1+(9)*((s)*(-2+(5)*(s))))*(-2+(2)*(s)+(5)*(t))));
  phit[24] = (2)*((sqrt(33))*((1+(9)*((s)*(-2+(5)*(s))))*(-3+(3)*(s)+(5)*(t)+(5)*(u))));
  phit[25] = (2)*((sqrt(55))*((1+(9)*((s)*(-2+(5)*(s))))*(-v+(3)*(u))));
  phit[26] = (10)*((sqrt(33))*((-1+(10)*(s))*(pow(-1+s,2)+(6)*((-1+s)*(t))+(7)*t2)));
  phit[27] = (2)*((sqrt(11))*((-1+(10)*(s))*(13+(13)*s2+(-24)*(u)+(s)*(-26+(66)*(t)+(24)*(u))+(3)*((t)*(-22+(21)*(t)+(28)*(u))))));
  phit[28] = (2)*((sqrt(165))*((-1+(10)*(s))*(3+(3)*s2+(7)*t2+(2)*((u)*(-8+(7)*(u)))+(2)*((s)*(-3+(5)*(t)+(8)*(u)))+(2)*((t)*(-5+(14)*(u))))));
  phit[29] = (2)*((sqrt(231))*((-1+(10)*(s))*(v2+(8)*((-v)*(u))+(10)*u2)));
  phit[30] = (4)*((sqrt(165))*((2)*(pow(-1+s,3))+(21)*((pow(-1+s,2))*(t))+(56)*((-1+s)*t2)+(42)*t3));
  phit[31] = (sqrt(55))*(-1+s3+(21)*(t)+(-84)*t2+(84)*t3+(3)*(s2*(-1+(7)*(t)))+(s)*(3+(42)*((t)*(-1+(2)*(t))))
             +(21)*((-1+s+(2)*(t))*((-1+s+(6)*(t))*(-v+(2)*(u)))));
  phit[32] = (5)*((sqrt(11/3.))*((2)*((pow(-1+s,2)+(16)*((-1+s)*(t))+(36)*t2)*(-v+(3)*(u)))+(8)*((-2+(2)*(s)+(9)*(t))*(v2+(6)*((-v)*(u))+(6)*u2))));
  phit[33] = (4)*((sqrt(1155))*(v2*(-1+s+(3)*(t))+(-v)*((-11+(11)*(s)+(27)*(t))*(u))+(5)*((-5+(5)*(s)+(9)*(t))*u2)+(15)*u3));
  phit[34] = (4)*((sqrt(165))*((-v3)+(15)*(v2*(u))+(45)*((-v)*u2)+(35)*u3));

  // phiu
  real_t* phiu = phi + 70;
  phiu[0] = 0;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = (2)*(sqrt(10));
  phiu[4] = 0;
  phiu[5] = 0;
  phiu[6] = (2)*((sqrt(14))*(-1+(6)*(s)));
  phiu[7] = 0;
  phiu[8] = (2)*((sqrt(21))*(-1+s+(5)*(t)));
  phiu[9] = (6)*((sqrt(35))*(-v+(2)*(u)));
  phiu[10] = 0;
  phiu[11] = 0;
  phiu[12] = (6)*((sqrt(2))*(1+(14)*((s)*(-1+(2)*(s)))));
  phiu[13] = 0;
  phiu[14] = (6)*((sqrt(3))*((-1+(8)*(s))*(-1+s+(5)*(t))));
  phiu[15] = (18)*((sqrt(5))*((-1+(8)*(s))*(-v+(2)*(u))));
  phiu[16] = 0;
  phiu[17] = (12)*(pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*(t2));
  phiu[18] = (12)*((sqrt(15))*((-1+s+(7)*(t))*(-v+(2)*(u))));
  phiu[19] = (24)*((sqrt(21))*(v2+(5)*((-v)*(u))+(5)*(u2)));
  phiu[20] = 0;
  phiu[21] = 0;
  phiu[22] = (2)*((sqrt(22))*(-1+(12)*((s)*(2+(s)*(-9+(10)*(s))))));
  phiu[23] = 0;
  phiu[24] = (2)*((sqrt(33))*((1+(9)*((s)*(-2+(5)*(s))))*(-1+s+(5)*(t))));
  phiu[25] = (6)*((sqrt(55))*((1+(9)*((s)*(-2+(5)*(s))))*(-v+(2)*(u))));
  phiu[26] = 0;
  phiu[27] = (4)*((sqrt(11))*((-1+(10)*(s))*(pow(-1+s,2)+(12)*((-1+s)*(t))+(21)*t2)));
  phiu[28] = (4)*((sqrt(165))*((-1+(10)*(s))*((-1+s+(7)*(t))*(-v+(2)*(u)))));
  phiu[29] = (8)*((sqrt(231))*((-1+(10)*(s))*(v2+(5)*((-v)*(u))+(5)*u2)));
  phiu[30] = 0;
  phiu[31] = (2)*((sqrt(55))*(pow(-1+s,3)+(21)*((pow(-1+s,2))*(t))+(84)*((-1+s)*t2)+(84)*t3));
  phiu[32] = (10)*((sqrt(33))*((pow(-1+s,2)+(16)*((-1+s)*(t))+(36)*t2)*(-v+(2)*(u))));
  phiu[33] = (4)*((sqrt(1155))*((-1+s+(9)*(t))*(v2+(5)*((-v)*(u))+(5)*u2)));
  phiu[34] = (20)*((sqrt(165))*((-v3)+(9)*(v2*(u))+(21)*((-v)*u2)+(14)*u3));
}

template<>
void
Legendre<Simplex,3,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * pentatopes, p = 5
*/
template<>
void
Legendre<Simplex,3,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Legendre<Simplex,3,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Legendre<Simplex,3,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro