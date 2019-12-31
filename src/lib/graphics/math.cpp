#include "common/error.h"

#include "graphics/math.h"

#include <numpack/dense/static/MatrixS_Transpose.h>
#include <numpack/Transpose.h>

namespace avro
{

namespace graphics
{

#if USE_GLM

#else

float*
mat3::value_ptr()
{
  int n = 0;
  for (int i=0;i<3;i++)
  for (int j=0;j<3;j++)
    mt_[n++] = (*this)(j,i);
  return mt_;
}

float*
mat4::value_ptr()
{
  int n = 0;
  for (int i=0;i<3;i++)
  for (int j=0;j<3;j++)
    mt_[n++] = (*this)(j,i);
  return mt_;
}

#endif

} // graphics

} // avro
