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

/*
 * pentatope, p = 1
*/
template<>
void
Bernstein<Simplex,4,1>::eval( const real_t* x , real_t* phi ) {
  real_t s = x[0];
  real_t t = x[1];
  real_t u = x[2];
  real_t v = x[3];

  phi[0] = 1 - s - t - u - v; // 1 at node 0 (s= t= u= v= 0)
  phi[1] =     s            ; // 1 at node 1 (s= 1, t= u= v= 0)
  phi[2] =         t        ; // 1 at node 2 (t= 1, s= u= v= 0)
  phi[3] =             u    ; // 1 at node 3 (u= 1, s= t= v= 0)
  phi[4] =                 v; // 1 at node 4 (v= 1, s= t= u= 0)
}

template<>
void
Bernstein<Simplex,4,1>::grad( const real_t* x , real_t* phi ) {

  real_t* phis = phi;
  phis[0] = -1;
  phis[1] = 1;
  phis[2] = 0;
  phis[3] = 0;
  phis[4] = 0;

  real_t* phit = phi + 5;
  phit[0] = -1;
  phit[1] = 0;
  phit[2] = 1;
  phit[3] = 0;
  phit[4] = 0;

  real_t* phiu = phi + 10;
  phiu[0] = -1;
  phiu[1] = 0;
  phiu[2] = 0;
  phiu[3] = 1;
  phiu[4] = 0;

  real_t *phiv = phi + 15;
  phiv[0] = -1;
  phiv[1] = 0;
  phiv[2] = 0;
  phiv[3] = 0;
  phiv[4] = 1;
}

template<>
void
Bernstein<Simplex,4,1>::hess( const real_t* x , real_t* phi ) {

  for (index_t i = 0; i < 50; i++)
    phi[i] = 0;
}

/*
 * pentatope, p = 2
*/
template<>
void
Bernstein<Simplex,4,2>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,2>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,2>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * pentatope, p = 3
*/
template<>
void
Bernstein<Simplex,4,3>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,3>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,3>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * pentatope, p = 4
*/
template<>
void
Bernstein<Simplex,4,4>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,4>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,4>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

/*
 * pentatope, p = 5
*/
template<>
void
Bernstein<Simplex,4,5>::eval( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,5>::grad( const real_t* x , real_t* phi ) {
  avro_implement;
}

template<>
void
Bernstein<Simplex,4,5>::hess( const real_t* x , real_t* phi ) {
  avro_implement;
}

} // avro
