#ifndef LUNA_LIB_NUMERICS_DETERMINANT_H_
#define LUNA_LIB_NUMERICS_DETERMINANT_H_

#include "common/types.h"

namespace luna
{

namespace numerics
{

template<typename type> class densMat;

template<typename type> type determinant(const densMat<type>& X);

} // numerics

} // luna

#endif
