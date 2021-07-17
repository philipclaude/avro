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

class vec3
to_vec3( const vec4& v) {
  vec3 x;
  x(0) = v(0);
  x(1) = v(1);
  x(2) = v(2);
  return x;
}

class vec4
to_vec4( const vec3& v , float h ) {
  vec4 x;
  x(0) = v(0);
  x(1) = v(1);
  x(2) = v(2);
  x(3) = h;
  return x;
}

class mat4
identity() {
  class mat4 m;
  for (index_t i = 0; i < 4; i++)
    m(i,i) = 1.0;
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
  float x = axis[0], y = axis[1], z = axis[2];
  float len = 1./norm(axis);

  x *= len;
  y *= len;
  z *= len;

  float s = sin(angle), c = cos(angle), t = 1 - c;

  float b00 = x * x * t + c;
  float b01 = y * x * t + z * s;
  float b02 = z * x * t - y * s;
  float b10 = x * y * t - z * s;
  float b11 = y * y * t + c;
  float b12 = z * y * t + x * s;
  float b20 = x * z * t + y * s;
  float b21 = y * z * t - x * s;
  float b22 = z * z * t + c;

  float a00 = a(0,0);
  float a01 = a(1,0);
  float a02 = a(2,0);
  float a03 = a(3,0);
  float a10 = a(0,1);
  float a11 = a(1,1);
  float a12 = a(2,1);
  float a13 = a(3,1);
  float a20 = a(0,2);
  float a21 = a(1,2);
  float a22 = a(2,2);
  float a23 = a(3,2);

  m(0,0) = a00 * b00 + a10 * b01 + a20 * b02;
  m(1,0) = a01 * b00 + a11 * b01 + a21 * b02;
  m(2,0) = a02 * b00 + a12 * b01 + a22 * b02;
  m(3,0) = a03 * b00 + a13 * b01 + a23 * b02;
  m(0,1) = a00 * b10 + a10 * b11 + a20 * b12;
  m(1,1) = a01 * b10 + a11 * b11 + a21 * b12;
  m(2,1) = a02 * b10 + a12 * b11 + a22 * b12;
  m(3,1) = a03 * b10 + a13 * b11 + a23 * b12;
  m(0,2) = a00 * b20 + a10 * b21 + a20 * b22;
  m(1,2) = a01 * b20 + a11 * b21 + a21 * b22;
  m(2,2) = a02 * b20 + a12 * b21 + a22 * b22;
  m(3,2) = a03 * b20 + a13 * b21 + a23 * b22;

  m(0,3) = a(0,3);
  m(1,3) = a(1,3);
  m(2,3) = a(2,3);
  m(3,3) = a(3,3);

  return m;
}

class mat4
translate( const class mat4& a , const vec3& t ) {
  class mat4 m;
  for (index_t i = 0; i < 4; i++)
  for (index_t j = 0; j < 4; j++)
    m(i,j) = a(i,j);

  m(0,3) = a(0,0)*t(0) + a(0,1)*t(1) + a(0,2)*t(2) + a(0,3);
  m(1,3) = a(1,0)*t(0) + a(1,1)*t(1) + a(1,2)*t(2) + a(1,3);
  m(2,3) = a(2,0)*t(0) + a(2,1)*t(1) + a(2,2)*t(2) + a(2,3);
  m(3,3) = a(3,0)*t(0) + a(3,1)*t(1) + a(3,2)*t(2) + a(3,3);

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
