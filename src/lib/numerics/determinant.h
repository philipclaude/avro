#ifndef luma_LIB_NUMERICS_DETERMINANT_H_
#define luma_LIB_NUMERICS_DETERMINANT_H_

#include "common/types.h"

namespace luma
{

namespace numerics
{

template<typename type> class densMat;

template<typename type> type determinant(const densMat<type>& X);

} // numerics

} // luma

#endif
