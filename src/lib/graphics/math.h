#ifndef URSA_LIB_GRAPHICS_MATH_H_
#define URSA_LIB_GRAPHICS_MATH_H_

#include "numerics/matrix.h"

namespace ursa
{

namespace graphics
{

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


} // graphics

} // ursa

#endif
