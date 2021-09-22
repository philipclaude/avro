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

#include <cmath>

namespace avro
{

/*
 * triangle, p = 0
*/
template<>
void
Legendre<Simplex,2,0>::eval( const real_t* x , real_t* phi ) {
  phi[0] = 1;
}

template<>
void
Legendre<Simplex,2,0>::grad( const real_t* x , real_t* phi ) {
  phi[0] = 0;
  phi[0] = 0;
}

template<>
void
Legendre<Simplex,2,0>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 3; i++)
    phi[i] = 0;
}

/*
 * triangle, p = 1
*/
template<>
void
Legendre<Simplex,2,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phi
  phi[0] = 1;
  phi[1] = (sqrt(2))*(-1+(3)*(s));
  phi[2] = (sqrt(6))*(-1+s+(2)*(t));
}

template<>
void
Legendre<Simplex,2,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (3)*(sqrt(2));
  phis[2] = sqrt(6);

  real_t* phit = phi + 3;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = (2)*(sqrt(6));
}

template<>
void
Legendre<Simplex,2,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 9; i++)
    phi[i] = 0;
}

/*
 * triangle, p = 2
*/
template<>
void
Legendre<Simplex,2,2>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t t2 = t*t;

  real_t u = s-1;
  real_t u2 = u*u;

  phi[0] = 1;
  phi[1] = (sqrt(2))*(-1+(3)*(s));
  phi[2] = (sqrt(6))*(-1+s+(2)*(t));
  phi[3] = (sqrt(3))*(1+(2)*((s)*(-4+(5)*(s))));
  phi[4] = (3)*((-1+(5)*(s))*(-1+s+(2)*(t)));
  phi[5] = (sqrt(15))*(u2+(6)*((-1+s)*(t))+(6)*(t2));
}

template<>
void
Legendre<Simplex,2,2>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (3)*(sqrt(2));
  phis[2] = sqrt(6);
  phis[3] = (4)*((sqrt(3))*(-2+(5)*(s)));
  phis[4] = (6)*(-3+(5)*(s)+(5)*(t));
  phis[5] = (2)*((sqrt(15))*(-1+s+(3)*(t)));

  real_t* phit = phi + 6;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = (2)*(sqrt(6));
  phit[3] = 0;
  phit[4] = -6+(30)*(s);
  phit[5] = (6)*((sqrt(15))*(-1+s+(2)*(t)));
}

template<>
void
Legendre<Simplex,2,2>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangle, p = 3
*/
template<>
void
Legendre<Simplex,2,3>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t t2 = t*t;

  real_t u = s-1;
  real_t u2 = u*u;

  phi[0] = 1;
  phi[1] = (sqrt(2))*(-1+(3)*(s));
  phi[2] = (sqrt(6))*(-1+s+(2)*(t));
  phi[3] = (sqrt(3))*(1+(2)*((s)*(-4+(5)*(s))));
  phi[4] = (3)*((-1+(5)*(s))*(-1+s+(2)*(t)));
  phi[5] = (sqrt(15))*(u2+(6)*((-1+s)*(t))+(6)*(t2));
  phi[6] = -2+(10)*((s)*(3+(s)*(-9+(7)*(s))));
  phi[7] = (2)*((sqrt(3))*((1+(3)*((s)*(-4+(7)*(s))))*(-1+s+(2)*(t))));
  phi[8] = (2)*((sqrt(5))*((-1+(7)*(s))*(u2+(6)*((-1+s)*(t))+(6)*(t2))));
  phi[9] = (2)*((sqrt(7))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
}

template<>
void
Legendre<Simplex,2,3>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t s2 = s*s;
  real_t t2 = t*t;
  real_t u = s-1;
  real_t u2 = u*u;

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (3)*(sqrt(2));
  phis[2] = sqrt(6);
  phis[3] = (4)*((sqrt(3))*(-2+(5)*(s)));
  phis[4] = (6)*(-3+(5)*(s)+(5)*(t));
  phis[5] = (2)*((sqrt(15))*(-1+s+(3)*(t)));
  phis[6] = (30)*(1+(s)*(-6+(7)*(s)));
  phis[7] = (2)*((sqrt(3))*(13+(-24)*(t)+(3)*((s)*(-22+(21)*(s)+(28)*(t)))));
  phis[8] = (6)*((sqrt(5))*(3+(7)*(s2)+(2)*((t)*(-8+(7)*(t)))+(2)*((s)*(-5+(14)*(t)))));
  phis[9] = (6)*((sqrt(7))*(u2+(8)*((-1+s)*(t))+(10)*(t2)));

  real_t* phit = phi + 10;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = (2)*(sqrt(6));
  phit[3] = 0;
  phit[4] = -6+(30)*(s);
  phit[5] = (6)*((sqrt(15))*(-1+s+(2)*(t)));
  phit[6] = 0;
  phit[7] = (4)*((sqrt(3))*(1+(3)*((s)*(-4+(7)*(s)))));
  phit[8] = (12)*((sqrt(5))*((-1+(7)*(s))*(-1+s+(2)*(t))));
  phit[9] = (24)*((sqrt(7))*(u2+(5)*((-1+s)*(t))+(5)*(t2)));
}

template<>
void
Legendre<Simplex,2,3>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangles, p = 4
*/
template<>
void
Legendre<Simplex,2,4>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  // phi
  real_t t2 = t*t;
  real_t t3 = t2*t;
  real_t t4 = t3*t;

  real_t u = s-1;
  real_t u2 = u*u;
  real_t u3 = u2*u;
  real_t u4 = u3*u;

  phi[0] = 1;
  phi[1] = (sqrt(2))*(-1+(3)*(s));
  phi[2] = (sqrt(6))*(-1+s+(2)*(t));
  phi[3] = (sqrt(3))*(1+(2)*((s)*(-4+(5)*(s))));
  phi[4] = (3)*((-1+(5)*(s))*(-1+s+(2)*(t)));
  phi[5] = (sqrt(15))*(u2+(6)*((-1+s)*(t))+(6)*(t2));
  phi[6] = -2+(10)*((s)*(3+(s)*(-9+(7)*(s))));
  phi[7] = (2)*((sqrt(3))*((1+(3)*((s)*(-4+(7)*(s))))*(-1+s+(2)*(t))));
  phi[8] = (2)*((sqrt(5))*((-1+(7)*(s))*(u2+(6)*((-1+s)*(t))+(6)*(t2))));
  phi[9] = (2)*((sqrt(7))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
  phi[10] = (sqrt(5))*(1+(2)*((s)*(-12+(7)*((s)*(9+(s)*(-16+(9)*(s)))))));
  phi[11] = (sqrt(15))*((-1+(21)*((pow(1+(-2)*(s),2))*(s)))*(-1+s+(2)*(t)));
  phi[12] = (5)*((1+(4)*((s)*(-4+(9)*(s))))*(u2+(6)*((-1+s)*(t))+(6)*(t2)));
  phi[13] = (sqrt(35))*((-1+(9)*(s))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
  phi[14] = (3)*((sqrt(5))*(u4+(20)*((u3)*(t))+(90)*((u2)*(t2))+(140)*((-1+s)*(t3))+(70)*(t4)));
}

template<>
void
Legendre<Simplex,2,4>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t s2 = s*s;
  real_t t2 = t*t;
  real_t t3 = t2*t;
  real_t u = s-1;
  real_t u2 = u*u;
  real_t u3 = u2*u;

  // phis
  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (3)*(sqrt(2));
  phis[2] = sqrt(6);
  phis[3] = (4)*((sqrt(3))*(-2+(5)*(s)));
  phis[4] = (6)*(-3+(5)*(s)+(5)*(t));
  phis[5] = (2)*((sqrt(15))*(-1+s+(3)*(t)));
  phis[6] = (30)*(1+(s)*(-6+(7)*(s)));
  phis[7] = (2)*((sqrt(3))*(13+(-24)*(t)+(3)*((s)*(-22+(21)*(s)+(28)*(t)))));
  phis[8] = (6)*((sqrt(5))*(3+(7)*(s2)+(2)*((t)*(-8+(7)*(t)))+(2)*((s)*(-5+(14)*(t)))));
  phis[9] = (6)*((sqrt(7))*(u2+(8)*((-1+s)*(t))+(10)*(t2)));
  phis[10] = (12)*((sqrt(5))*(-2+(7)*((s)*(3+(-8)*(s)+(6)*(s2)))));
  phis[11] = (2)*((sqrt(15))*(-11+(21)*(t)+(21)*((s)*(5+(-8)*(t)+(4)*((s)*(-3+(2)*(s)+(3)*(t)))))));
  phis[12] = (30)*(-3+(s)*(23+(4)*((s)*(-11+(6)*(s))))+(17)*(t)+(4)*((s)*((-26+(27)*(s))*(t)))+(8)*((-2+(9)*(s))*(t2)));
  phis[13] = (12)*((sqrt(35))*((u2)*(-1+(3)*(s))+(-1+s)*((-11+(27)*(s))*(t))+(5)*((-5+(9)*(s))*(t2))+(15)*(t3)));
  phis[14] = (12)*((sqrt(5))*(u3+(15)*((u2)*(t))+(45)*((-1+s)*(t2))+(35)*(t3)));

  // phit
  real_t *phit = phi + 15;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = (2)*(sqrt(6));
  phit[3] = 0;
  phit[4] = -6+(30)*(s);
  phit[5] = (6)*((sqrt(15))*(-1+s+(2)*(t)));
  phit[6] = 0;
  phit[7] = (4)*((sqrt(3))*(1+(3)*((s)*(-4+(7)*(s)))));
  phit[8] = (12)*((sqrt(5))*((-1+(7)*(s))*(-1+s+(2)*(t))));
  phit[9] = (24)*((sqrt(7))*(u2+(5)*((-1+s)*(t))+(5)*(t2)));
  phit[10] = 0;
  phit[11] = (2)*((sqrt(15))*(-1+(21)*((pow(1+(-2)*(s),2))*(s))));
  phit[12] = (30)*((1+(4)*((s)*(-4+(9)*(s))))*(-1+s+(2)*(t)));
  phit[13] = (12)*((sqrt(35))*((-1+(9)*(s))*(u2+(5)*((-1+s)*(t))+(5)*(t2))));
  phit[14] = (60)*((sqrt(5))*((-1+s+(2)*(t))*(u2+(7)*((-1+s)*(t))+(7)*(t2))));
}

template<>
void
Legendre<Simplex,2,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * triangles, p = 5
*/
template<>
void
Legendre<Simplex,2,5>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t s2 = s*s;

  real_t t2 = t*t;
  real_t t3 = t2*t;
  real_t t4 = t3*t;
  real_t t5 = t4*t;

  real_t u = s-1;
  real_t u2 = u*u;
  real_t u3 = u2*u;
  real_t u4 = u3*u;
  real_t u5 = u4*u;

  phi[0] = 1;
  phi[1] = (sqrt(2))*(-1+(3)*(s));
  phi[2] = (sqrt(6))*(-1+s+(2)*(t));
  phi[3] = (sqrt(3))*(1+(2)*((s)*(-4+(5)*(s))));
  phi[4] = (3)*((-1+(5)*(s))*(-1+s+(2)*(t)));
  phi[5] = (sqrt(15))*(u2+(6)*((-1+s)*(t))+(6)*(t2));
  phi[6] = -2+(10)*((s)*(3+(s)*(-9+(7)*(s))));
  phi[7] = (2)*((sqrt(3))*((1+(3)*((s)*(-4+(7)*(s))))*(-1+s+(2)*(t))));
  phi[8] = (2)*((sqrt(5))*((-1+(7)*(s))*(u2+(6)*((-1+s)*(t))+(6)*(t2))));
  phi[9] = (2)*((sqrt(7))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
  phi[10] = (sqrt(5))*(1+(2)*((s)*(-12+(7)*((s)*(9+(s)*(-16+(9)*(s)))))));
  phi[11] = (sqrt(15))*((-1+(21)*((pow(1+(-2)*(s),2))*(s)))*(-1+s+(2)*(t)));
  phi[12] = (5)*((1+(4)*((s)*(-4+(9)*(s))))*(u2+(6)*((-1+s)*(t))+(6)*(t2)));
  phi[13] = (sqrt(35))*((-1+(9)*(s))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
  phi[14] = (3)*((sqrt(5))*(u4+(20)*((u3)*(t))+(90)*((u2)*(t2))+(140)*((-1+s)*(t3))+(70)*(t4)));
  phi[15] = (sqrt(6))*(-1+(7)*((s)*(5+(2)*((s)*(-20+(3)*((s)*(20+(s)*(-25+(11)*(s)))))))));
  phi[16] = (3)*((sqrt(2))*((1+(2)*((s)*(-16+(3)*((s)*(36+(-80)*(s)+(55)*(s2))))))*(-1+s+(2)*(t))));
  phi[17] = (sqrt(30))*((-1+(3)*((s)*(9+(5)*((s)*(-9+(11)*(s))))))*(u2+(6)*((-1+s)*(t))+(6)*(t2)));
  phi[18] = (sqrt(42))*((1+(5)*((s)*(-4+(11)*(s))))*((-1+s+(2)*(t))*(u2+(10)*((-1+s)*(t))+(10)*(t2))));
  phi[19] = (3)*((sqrt(6))*((-1+(11)*(s))*(u4+(20)*((u3)*(t))+(90)*((u2)*(t2))+(140)*((-1+s)*(t3))
            +(70)*(t4))));
  phi[20] = (sqrt(66))*(u5+(30)*((u4)*(t))+(210)*((u3)*(t2))+(560)*((u2)*(t3))
            +(630)*((-1+s)*(t4))+(252)*(t5));
}

template<>
void
Legendre<Simplex,2,5>::grad( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];

  real_t s2 = s*s;

  real_t t2 = t*t;
  real_t t3 = t2*t;
  real_t t4 = t3*t;

  real_t u = s-1;
  real_t u2 = u*u;
  real_t u3 = u2*u;
  real_t u4 = u3*u;

  real_t* phis = phi;
  phis[0] = 0;
  phis[1] = (3)*(sqrt(2));
  phis[2] = sqrt(6);
  phis[3] = (4)*((sqrt(3))*(-2+(5)*(s)));
  phis[4] = (6)*(-3+(5)*(s)+(5)*(t));
  phis[5] = (2)*((sqrt(15))*(-1+s+(3)*(t)));
  phis[6] = (30)*(1+(s)*(-6+(7)*(s)));
  phis[7] = (2)*((sqrt(3))*(13+(-24)*(t)+(3)*((s)*(-22+(21)*(s)+(28)*(t)))));
  phis[8] = (6)*((sqrt(5))*(3+(7)*(s2)+(2)*((t)*(-8+(7)*(t)))+(2)*((s)*(-5+(14)*(t)))));
  phis[9] = (6)*((sqrt(7))*(u2+(8)*((-1+s)*(t))+(10)*(t2)));
  phis[10] = (12)*((sqrt(5))*(-2+(7)*((s)*(3+(-8)*(s)+(6)*(s2)))));
  phis[11] = (2)*((sqrt(15))*(-11+(21)*(t)+(21)*((s)*(5+(-8)*(t)+(4)*((s)*(-3+(2)*(s)+(3)*(t)))))));
  phis[12] = (30)*(-3+(s)*(23+(4)*((s)*(-11+(6)*(s))))+(17)*(t)+(4)*((s)*((-26+(27)*(s))*(t)))+(8)*((-2+(9)*(s))*(t2)));
  phis[13] = (12)*((sqrt(35))*((u2)*(-1+(3)*(s))+(-1+s)*((-11+(27)*(s))*(t))+(5)*((-5+(9)*(s))*(t2))+(15)*(t3)));
  phis[14] = (12)*((sqrt(5))*(u3+(15)*((u2)*(t))+(45)*((-1+s)*(t2))+(35)*(t3)));
  phis[15] = (35)*((sqrt(6))*(1+(2)*((s)*(-8+(3)*((s)*(12+(s)*(-20+(11)*(s))))))));
  phis[16] = (3)*((sqrt(2))*(33+(-64)*(t)+(2)*((s)*((8)*(-31+(54)*(t))+(3)*((s)*(348+(-480)*(t)+(5)*((s)*(-108+(55)*(s)+(88)*(t)))))))));
  phis[17] = (sqrt(30))*((2)*((-1+(3)*((s)*(9+(5)*((s)*(-9+(11)*(s))))))*(-1+s+(3)*(t)))+(9)*((3+(5)*((s)*(-6+(11)*(s))))*(u2
             +(6)*((-1+s)*(t))+(6)*(t2))));
  phis[18] = (sqrt(42))*((u2)*(23+(5)*((s)*(-38+(55)*(s))))+(24)*((-1+s)*((11+(5)*((s)*(-17+(22)*(s))))*(t)))
             +(90)*((7+(-50)*(s)+(55)*(s2))*(t2))+(200)*((-2+(11)*(s))*(t3)));
  phis[19] = (15)*((sqrt(6))*((u3)*(-3+(11)*(s))+(8)*((u2)*((-7+(22)*(s))*(t)))+(18)*((-1+s)*((-13+(33)*(s))*(t2)))
             +(56)*((-6+(11)*(s))*(t3))+(154)*(t4)));
  phis[20] = (5)*((sqrt(66))*(u4+(24)*((u3)*(t))+(126)*((u2)*(t2))+(224)*((-1+s)*(t3))+(126)*(t4)));

  real_t *phit = phi + 21;
  phit[0] = 0;
  phit[1] = 0;
  phit[2] = (2)*(sqrt(6));
  phit[3] = 0;
  phit[4] = -6+(30)*(s);
  phit[5] = (6)*((sqrt(15))*(-1+s+(2)*(t)));
  phit[6] = 0;
  phit[7] = (4)*((sqrt(3))*(1+(3)*((s)*(-4+(7)*(s)))));
  phit[8] = (12)*((sqrt(5))*((-1+(7)*(s))*(-1+s+(2)*(t))));
  phit[9] = (24)*((sqrt(7))*(u2+(5)*((-1+s)*(t))+(5)*(t2)));
  phit[10] = 0;
  phit[11] = (2)*((sqrt(15))*(-1+(21)*((pow(1+(-2)*(s),2))*(s))));
  phit[12] = (30)*((1+(4)*((s)*(-4+(9)*(s))))*(-1+s+(2)*(t)));
  phit[13] = (12)*((sqrt(35))*((-1+(9)*(s))*(u2+(5)*((-1+s)*(t))+(5)*(t2))));
  phit[14] = (60)*((sqrt(5))*((-1+s+(2)*(t))*(u2+(7)*((-1+s)*(t))+(7)*(t2))));
  phit[15] = 0;
  phit[16] = (6)*((sqrt(2))*(1+(2)*((s)*(-16+(3)*((s)*(36+(-80)*(s)+(55)*(s2)))))));
  phit[17] = (6)*((sqrt(30))*((-1+(3)*((s)*(9+(5)*((s)*(-9+(11)*(s))))))*(-1+s+(2)*(t))));
  phit[18] = (12)*((sqrt(42))*((1+(5)*((s)*(-4+(11)*(s))))*(u2+(5)*((-1+s)*(t))+(5)*(t2))));
  phit[19] = (60)*((sqrt(6))*((-1+(11)*(s))*((-1+s+(2)*(t))*(u2+(7)*((-1+s)*(t))+(7)*(t2)))));
  phit[20] = (30)*((sqrt(66))*(u4+(14)*((u3)*(t))+(56)*((u2)*(t2))+(84)*((-1+s)*(t3))+(42)*(t4)));
}

template<>
void
Legendre<Simplex,2,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
