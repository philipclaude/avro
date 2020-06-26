//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_MATH_H_
#define avro_LIB_GRAPHICS_MATH_H_

#include "numerics/matrix.h"

#define USE_GLM 1

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

namespace glm
{
typedef numerics::VectorS<2,float> vec2;
typedef numerics::VectorS<3,float> vec3;
typedef numerics::VectorS<4,float> vec4;

class mat3 : public numerics::MatrixS<3,3,float>
{
public:
  float* value_ptr();

  float* operator[] (index_t i) { return &(*this)(i,0); }
};

class mat4 : public numerics::MatrixS<4,4,float>
{
public:
  float* value_ptr();
  float* operator[] (index_t i) { return &(*this)(i,0); }
};

} // glm

using mat4 = glm::mat4;

#endif


} // graphics

} // avro

#endif
