//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_MATH_H_
#define avro_LIB_GRAPHICS_MATH_H_

#include "numerics/mat.h"
#include "numerics/mat.hpp"
#include "numerics/vec.h"

#define USE_GLM 0

#if USE_GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

#endif

namespace avro
{

namespace graphics
{


#if USE_GLM

using vec2 = glm::vec2;
using vec3 = glm::vec3;
using vec4 = glm::vec4;

using mat3 = glm::mat3;
using mat4 = glm::mat4;

#else

class vec2 : public vecs<2,float> {
public:
  vec2() {}

  vec2( const vecs<2,float>& v ) {
    for (index_t i = 0; i < 2; i++)
      (*this)(i) = v(i);
  }

  template<typename R>
  vec2( const std::initializer_list<R>& v ) {
    float v0 = float( *(v.begin()+0) );
    float v1 = float( *(v.begin()+1) );
    operator=( {v0,v1} );
  }

  vec2& operator= (const std::initializer_list<float>& v ) {
    avro_assert( v.size() == 2 );
    index_t i = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      (*this)(i++) = *it;
    return *this;
  }
};

class vec3 : public vecs<3,float> {
public:
  vec3() {}

  vec3( const vecs<3,float>& v ) {
    for (index_t i = 0; i < 3; i++)
      (*this)(i) = v(i);
  }

  template<typename R>
  vec3( const std::initializer_list<R>& v ) {
    float v0 = float( *(v.begin()+0) );
    float v1 = float( *(v.begin()+1) );
    float v2 = float( *(v.begin()+2) );
    operator=( {v0,v1,v2} );
  }

  template<typename R>
  vec3( const std::vector<R>& v ) {
    avro_assert( v.size() == 3 );
    (*this)(0) = v[0];
    (*this)(1) = v[1];
    (*this)(2) = v[2];
  }

  vec3& operator= (const std::initializer_list<float>& v ) {
    avro_assert( v.size() == 3 );
    index_t i = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      (*this)(i++) = *it;
    return *this;
  }
};

class vec4 : public vecs<4,float> {
public:
  vec4() : vecs<4,float>() {}

  vec4( const vecs<4,float>& v ) {
    for (index_t i = 0; i < 4; i++)
      (*this)(i) = v(i);
  }

  template<typename R>
  vec4( const std::initializer_list<R>& v ) {
    float v0 = float( *(v.begin()+0) );
    float v1 = float( *(v.begin()+1) );
    float v2 = float( *(v.begin()+2) );
    float v3 = float( *(v.begin()+3) );
    operator=( {v0,v1,v2,v3} );
  }

  vec4& operator= (const std::initializer_list<float>& v ) {
    avro_assert( v.size() == 4 );
    index_t i = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      (*this)(i++) = *it;
    return *this;
  }
};

class mat3 : public mats<3,3,float> {
public:

  mat3() :
    mats<3,3,float>()
  {}

  mat3( const mats<3,3,float>& m ) {
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*this)(i,j) = m(i,j);
  }

  float* operator[] (index_t i) { return &(*this)(i,0); }
};

class mat4 : public mats<4,4,float> {
public:
  mat4() :
    mats<4,4,float>()
  {}

  mat4( const mats<4,4,float>& m ) {
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      (*this)(i,j) = m(i,j);
  }

  template<typename R>
  mat4( const std::vector<R>& m ) {
    avro_assert( m.size() == 16 );
    for (int i = 0; i < 16; i++)
      data_[i] = m[i];
  }

  float* operator[] (index_t i) { return &(*this)(i,0); }
  const float* operator[] (index_t i) const { return &(*this)(i,0); }
};

namespace glm {

class mat4 mat4( float a );
class vec3 to_vec3( const vec4& x );
class vec4 to_vec4( const vec3& x , float h );
class mat4 identity();
class mat4 perspective( float fov , float aspect , float znear , float zfar );
class mat4 lookAt( const vec3& eye , const vec3& center , const vec3& up );
class mat4 rotate( const class mat4& m , float angle , const vec3& axis );
class mat4 translate( const class mat4& m , const vec3& t );
class mat4 scale( const class mat4& m , const vec3& s );
class mat4 inverse( const class mat4& m );
class mat3 inverse( const class mat3& m );
class mat4 transpose( const class mat4& m );
std::string to_string( const class mat4& m );
class vec3 cross( const vec3& u , const vec3& v );
float norm( const vec3& u );
class vec3 normalize( const vec3& u );
float dot( const vec3& u , const vec3& v);

}

// component-wise vector multiplication (use in ray tracing)
vec3 operator* (const vec3& u , const vec3& v);

#endif


} // graphics

} // avro

#endif
