#ifndef LUNA_LIB_GRAPHICS_MATH_H_
#define LUNA_LIB_GRAPHICS_MATH_H_

#include "numerics/matrix.h"

#define USE_GLM 1

#if USE_GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#endif

namespace luna
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

typedef numerics::VectorS<2,float> vec2;
typedef numerics::VectorS<3,float> vec3;
typedef numerics::VectorS<4,float> vec4;

class mat3 : public numerics::MatrixS<3,3,float>
{
public:
  float* value_ptr();

private:
  float mt_[9];
};

class mat4 : public numerics::MatrixS<4,4,float>
{
public:
  float* value_ptr();
private:
  float mt_[16];
};

#endif


} // graphics

} // luna

#endif
