//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"

#include "graphics/math.h"

#include "numerics/linear_algebra.h"
#include "numerics/vec.hpp"

namespace avro
{

namespace graphics
{

#if USE_GLM == 0

namespace glm
{

class mat4
mat4( float a ) {
  class mat4 m;
  for (index_t i = 0; i < 4; i++)
    m(i,i) = a;
  return m;
}

class mat4
perspective( float fov , float aspect , float n , float f ) {
  class mat4 m;

  float a = 1./tan(fov/2.0);

  m(0,0) = a / aspect;
  m(1,1) = a;
  m(2,2) = (f+n)/(n-f);
  m(2,3) = 2*f*n/(n-f);
  m(3,2) = -1.0;

  return m;
}

vec3
cross( const vec3& u , const vec3& v ) {
  vec3 w;
  w(0) =    u(1)*v(2) - u(2)*v(1);
  w(1) = -( u(0)*v(2) - u(2)*v(0) );
  w(2) =    u(0)*v(1) - u(1)*v(0);
  return w;
}

float
norm( const vec3& u ) {
  return std::sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] );
}

class mat4
lookAt( const vec3& eye , const vec3& center , const vec3& up ) {
  class mat4 m;

  // calculate the gaze
  vec3 g = center - eye;
  vec3 w = g * real_t(-1./norm(g));
  vec3 u = cross(up,w);
  u = u * (1./norm(u));
  vec3 v = cross(w,u);

  m(0,0) = u(0);
  m(1,0) = v(0);
  m(2,0) = w(0);
  m(0,1) = u(1);
  m(1,1) = v(1);
  m(2,1) = w(1);
  m(0,2) = u(2);
  m(1,2) = v(2);
  m(2,2) = w(2);

  m(0,3) = -eye(0)*u(0) - eye(1)*u(1) - eye(2)*u(2);
  m(1,3) = -eye(0)*v(0) - eye(1)*v(1) - eye(2)*v(2);
  m(2,3) = -eye(0)*w(0) - eye(1)*w(1) - eye(2)*w(2);
  m(3,3) = 1.0;
  return m;
}

class mat4
scale( const class mat4& a , const vec3& s ) {
  class mat4 m;
  for (index_t i = 0; i < 4; i++)
  for (index_t j = 0; j < 3; j++)
    m(i,j) = a(i,j)*s(j);

  for (index_t i = 0; i < 4; i++)
    m(i,3) = a(i,3);
  return m;
}

class mat4
rotate( const class mat4& a , float angle , const vec3& axis ) {
  class mat4 m;
  avro_implement;
  return m;
}

class mat4
translate( const class mat4& a , const vec3& t ) {
  class mat4 m;
  avro_implement;
  return m;
}

class mat4
transpose( const class mat4& a ) {
  return numerics::transpose(a);
}

class mat4
inverse( const class mat4& a ) {
  return numerics::inverse(a);
}

std::string
to_string( const class mat4& m ) {
  std::string s = "{ ";
  for (index_t i = 0; i < 4; i++)
  for (index_t j = 0; j < 4; j++)
    s += std::to_string(m(i,j)) + " ";
  s += "}";
  return s;
}

} // glm


#endif

} // graphics

} // avro
